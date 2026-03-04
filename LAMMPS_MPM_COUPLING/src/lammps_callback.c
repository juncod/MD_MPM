/*
 * lammps_callback.c — LAMMPS fix external callback for MPM coupling
 *
 * Design:
 *   - All physics parameters defined as LAMMPS variables in .in file
 *   - C driver reads them via lammps_extract_variable()
 *   - Callback frequency controlled by fix external Nevery (not C code)
 *   - run N controlled by .in file variable, C issues the command
 *
 * Callback flow (invoked by LAMMPS every Nevery steps):
 *   1. MPI_Gatherv: all rank atom positions + CNA → rank 0
 *   2. Rank 0: atom→MP, MPM solve -div(σ∇φ)=0, compute F=q*E
 *   3. MPI_Scatterv: forces back to all ranks
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "lammps_callback.h"
#include "atom_to_mp.h"
#include "mpm_electric.h"
#include "mpm_mesh.h"
#include "mpm_solver.h"
#include "vtk_output.h"

/* LAMMPS library interface */
#include "library.h"

/* --------------------------------------------------------------------------
 * Read a LAMMPS equal-style variable, return its value (0.0 if not found)
 * --------------------------------------------------------------------------*/
static double read_lammps_variable(void *lmp, const char *name, double fallback)
{
    char *val = (char *)lammps_extract_variable(lmp, (char *)name, NULL);
    if (val) {
        double result = atof(val);
        /* LAMMPS allocates the string; must free via lammps_free */
        lammps_free(val);
        return result;
    }
    fprintf(stderr, "Warning: LAMMPS variable '%s' not found, using %.6g\n",
            name, fallback);
    return fallback;
}

/* --------------------------------------------------------------------------
 * Force callback — called by LAMMPS every Nevery steps
 *
 * Signature for fix external pf/callback:
 *   void callback(void *caller, int64_t timestep, int nlocal,
 *                 int *ids, double **x, double **f)
 *
 * Note: Nevery is set in the .in file:
 *   fix mpm_efield mobile external pf/callback 100 3
 *   → LAMMPS calls this callback every 100 steps automatically
 * --------------------------------------------------------------------------*/
void mpm_force_callback(void *caller, int64_t timestep,
                        int nlocal, int *ids, double **x, double **f)
{
    CouplingState *state = (CouplingState *)caller;

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* --- Gather atom data to rank 0 --- */
    int *recv_counts = NULL;
    int *displs = NULL;
    int natoms_total = 0;

    if (rank == 0) {
        recv_counts = (int *)mpm_calloc(nprocs, sizeof(int));
        displs = (int *)mpm_calloc(nprocs, sizeof(int));
    }
    MPI_Gather(&nlocal, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < nprocs; i++) {
            displs[i] = natoms_total;
            natoms_total += recv_counts[i];
        }
    }
    MPI_Bcast(&natoms_total, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Flatten local positions */
    double *local_pos = (double *)mpm_calloc(nlocal * 3, sizeof(double));
    for (int i = 0; i < nlocal; i++) {
        local_pos[3 * i + 0] = x[i][0];
        local_pos[3 * i + 1] = x[i][1];
        local_pos[3 * i + 2] = x[i][2];
    }

    /* Gather positions (×3) and IDs */
    int *recv_counts3 = NULL, *displs3 = NULL;
    double *all_pos = NULL;
    int *all_ids = NULL;

    if (rank == 0) {
        recv_counts3 = (int *)mpm_calloc(nprocs, sizeof(int));
        displs3 = (int *)mpm_calloc(nprocs, sizeof(int));
        for (int i = 0; i < nprocs; i++) {
            recv_counts3[i] = recv_counts[i] * 3;
            displs3[i] = displs[i] * 3;
        }
        all_pos = (double *)mpm_calloc(natoms_total * 3, sizeof(double));
        all_ids = (int *)mpm_calloc(natoms_total, sizeof(int));
    }

    MPI_Gatherv(local_pos, nlocal * 3, MPI_DOUBLE,
                 all_pos, recv_counts3, displs3, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    MPI_Gatherv(ids, nlocal, MPI_INT,
                 all_ids, recv_counts, displs, MPI_INT,
                 0, MPI_COMM_WORLD);

    /* Gather CNA values */
    double *cna_ptr = (double *)lammps_extract_compute(
        state->lmp, "cna_mpm", 1, 1);

    int *local_cna = (int *)mpm_calloc(nlocal, sizeof(int));
    if (cna_ptr) {
        for (int i = 0; i < nlocal; i++)
            local_cna[i] = (int)cna_ptr[i];
    } else {
        for (int i = 0; i < nlocal; i++)
            local_cna[i] = 1;  /* assume FCC if CNA not available */
    }

    int *all_cna = NULL;
    if (rank == 0)
        all_cna = (int *)mpm_calloc(natoms_total, sizeof(int));

    MPI_Gatherv(local_cna, nlocal, MPI_INT,
                 all_cna, recv_counts, displs, MPI_INT,
                 0, MPI_COMM_WORLD);

    /* --- Rank 0: MPM solve --- */
    double *all_forces = NULL;
    if (rank == 0) {
        all_forces = (double *)mpm_calloc(natoms_total * 3, sizeof(double));

        /* Convert atoms to material points */
        atoms_to_mps(&state->mpdata, &state->config, all_pos, all_cna,
                     all_ids, natoms_total);

        /* Reset solver for fresh solve */
        memset(state->solver.phi, 0, state->solver.n * sizeof(double));
        memset(state->solver.frct, 0, state->solver.n * sizeof(double));

        /* Solve electric field */
        int ret = mpm_electric_solve(state);
        if (ret != 0) {
            fprintf(stderr, "MPM solve failed at step %lld\n",
                    (long long)timestep);
        }

        /* Extract forces */
        mps_to_forces(&state->mpdata, all_forces, natoms_total);

        if (timestep % 10000 == 0) {
            fprintf(stdout, "  [MPM] step %lld: %d atoms, solved\n",
                    (long long)timestep, natoms_total);
        }
    }

    /* --- Scatter forces back --- */
    double *local_forces = (double *)mpm_calloc(nlocal * 3, sizeof(double));
    MPI_Scatterv(all_forces, recv_counts3, displs3, MPI_DOUBLE,
                 local_forces, nlocal * 3, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    /* Apply forces to fix external array */
    for (int i = 0; i < nlocal; i++) {
        f[i][0] += local_forces[3 * i + 0];
        f[i][1] += local_forces[3 * i + 1];
        f[i][2] += local_forces[3 * i + 2];
    }

    /* Cleanup */
    free(local_pos);
    free(local_cna);
    free(local_forces);
    if (rank == 0) {
        free(recv_counts);  free(displs);
        free(recv_counts3); free(displs3);
        free(all_pos);      free(all_ids);
        free(all_cna);      free(all_forces);
    }
}

/* --------------------------------------------------------------------------
 * Initialize: open LAMMPS, execute .in (up to run 0), read variables,
 *             build mesh, register callback
 * --------------------------------------------------------------------------*/
int coupling_init(CouplingState *state, const char *input_script)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Execute input script (everything including "run 0") */
    lammps_file(state->lmp, (char *)input_script);

    /* Read coupling parameters from LAMMPS variables */
    MPMConfig *cfg = &state->config;
    cfg->voltage_left  = read_lammps_variable(state->lmp, "voltage_left",  0.0);
    cfg->voltage_right = read_lammps_variable(state->lmp, "voltage_right", 1.0);
    cfg->sigma_bulk    = read_lammps_variable(state->lmp, "sigma_bulk",    1.0);
    cfg->sigma_defect  = read_lammps_variable(state->lmp, "sigma_defect",  0.1);
    cfg->atom_charge   = read_lammps_variable(state->lmp, "atom_charge",   1.0);
    cfg->target_h      = read_lammps_variable(state->lmp, "target_h",     10.0);
    cfg->nr_max_iter   = 5;
    cfg->nr_tolerance  = 1e-9;
    cfg->mp_type       = 2;  /* GIMP */

    /* Extract simulation box */
    double boxlo[3], boxhi[3], xy, yz, xz;
    int pflags[3], boxflag;
    lammps_extract_box(state->lmp, boxlo, boxhi, &xy, &yz, &xz, pflags, &boxflag);

    cfg->x0 = boxlo[0];  cfg->y0 = boxlo[1];  cfg->z0 = boxlo[2];
    cfg->lx = boxhi[0] - boxlo[0];
    cfg->ly = boxhi[1] - boxlo[1];
    cfg->lz = boxhi[2] - boxlo[2];

    /* Auto-compute mesh from target_h */
    double th = cfg->target_h;
    if (th <= 0.0) th = 10.0;
    cfg->nelsx = (int)round(cfg->lx / th);
    cfg->nelsy = (int)round(cfg->ly / th);
    cfg->nelsz = (int)round(cfg->lz / th);
    if (cfg->nelsx < 1) cfg->nelsx = 1;
    if (cfg->nelsy < 1) cfg->nelsy = 1;
    if (cfg->nelsz < 1) cfg->nelsz = 1;

    if (rank == 0) {
        printf("\n=== MPM Coupling Parameters (from .in) ===\n");
        printf("  Voltage: left=%.4f, right=%.4f\n", cfg->voltage_left, cfg->voltage_right);
        printf("  Conductivity: bulk=%.4f, defect=%.4f\n", cfg->sigma_bulk, cfg->sigma_defect);
        printf("  Atom charge: %.6e\n", cfg->atom_charge);
        printf("  Box: [%.1f, %.1f, %.1f] to [%.1f, %.1f, %.1f]\n",
               boxlo[0], boxlo[1], boxlo[2], boxhi[0], boxhi[1], boxhi[2]);
        printf("  Domain: %.1f x %.1f x %.1f Angstrom\n", cfg->lx, cfg->ly, cfg->lz);
        printf("  Mesh: %d x %d x %d elements (h ~ %.1f x %.1f x %.1f)\n",
               cfg->nelsx, cfg->nelsy, cfg->nelsz,
               cfg->lx / cfg->nelsx, cfg->ly / cfg->nelsy, cfg->lz / cfg->nelsz);
    }

    /* Create mesh and solver (rank 0 only) */
    if (rank == 0) {
        mesh_create(&state->mesh, cfg);
        mesh_set_dirichlet_bc(&state->mesh, cfg);
        solver_init(&state->solver, state->mesh.nnodes);
        printf("  MPM mesh: %d nodes, %d elements\n",
               state->mesh.nnodes, state->mesh.nels);
        printf("==========================================\n\n");
    }

    /* Register callback on fix "mpm_efield" */
    lammps_set_fix_external_callback(state->lmp, "mpm_efield",
                                     mpm_force_callback, state);

    state->initialized = 1;
    return 0;
}

/* --------------------------------------------------------------------------
 * Run: read run_steps from LAMMPS variable, issue "run N"
 * --------------------------------------------------------------------------*/
int coupling_run(CouplingState *state)
{
    int nsteps = (int)read_lammps_variable(state->lmp, "run_steps", 10000);
    int vtk_freq = (int)read_lammps_variable(state->lmp, "vtk_freq", 0);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (vtk_freq > 0) {
        /* Run in chunks for VTK output */
        int steps_done = 0;
        int vtk_step = 0;
        while (steps_done < nsteps) {
            int chunk = vtk_freq;
            if (steps_done + chunk > nsteps)
                chunk = nsteps - steps_done;

            char cmd[256];
            snprintf(cmd, sizeof(cmd), "run %d", chunk);
            lammps_command(state->lmp, cmd);

            steps_done += chunk;
            vtk_step++;

            if (rank == 0 && state->mpdata.nmp > 0) {
                vtk_write_step(state, "output", vtk_step);
            }
        }
    } else {
        /* Single run, no VTK */
        char cmd[256];
        snprintf(cmd, sizeof(cmd), "run %d", nsteps);
        lammps_command(state->lmp, cmd);
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Cleanup
 * --------------------------------------------------------------------------*/
void coupling_free(CouplingState *state)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        mpdata_free(&state->mpdata);
        solver_free(&state->solver);
        mesh_free(&state->mesh);
    }

    free(state->prev_forces);
    state->prev_forces = NULL;

    if (state->lmp) {
        lammps_close(state->lmp);
        state->lmp = NULL;
    }
}
