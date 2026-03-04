/*
 * main.c — LAMMPS-MPM Electric Field Coupling Driver
 *
 * Couples LAMMPS molecular dynamics with MPM electric field solver.
 * At every coupling_interval steps, solves -div(sigma*grad(phi))=0
 * using atom positions as material points, then applies F=q*E forces.
 *
 * Usage: mpirun -np N ./lammps_mpm_coupling [options] input.in
 *
 * Options:
 *   -nelsx N      — elements in x (default: 20)
 *   -nelsy N      — elements in y (default: 20)
 *   -nelsz N      — elements in z (default: 20)
 *   -sigma_bulk V — bulk conductivity (default: 1.0)
 *   -sigma_def V  — defect conductivity (default: 0.1)
 *   -vleft V      — left BC voltage (default: 0.0)
 *   -vright V     — right BC voltage (default: 1.0)
 *   -charge Q     — atom charge for F=q*E (default: 1.0)
 *   -interval N   — coupling interval in steps (default: 100)
 *   -nsteps N     — total MD steps (default: 10000)
 *   -vtk_dir DIR  — VTK output directory (default: output)
 *   -vtk_freq N   — VTK output frequency (default: 1000)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "mpm_types.h"
#include "lammps_callback.h"
#include "vtk_output.h"

static void print_usage(const char *prog)
{
    fprintf(stderr,
        "Usage: mpirun -np N %s [options] input.in\n"
        "Options:\n"
        "  -nelsx N      elements in x (default: 20)\n"
        "  -nelsy N      elements in y (default: 20)\n"
        "  -nelsz N      elements in z (default: 20)\n"
        "  -sigma_bulk V bulk conductivity (default: 1.0)\n"
        "  -sigma_def V  defect conductivity (default: 0.1)\n"
        "  -vleft V      left BC voltage (default: 0.0)\n"
        "  -vright V     right BC voltage (default: 1.0)\n"
        "  -charge Q     atom charge for F=q*E (default: 1.0)\n"
        "  -interval N   coupling interval (default: 100)\n"
        "  -nsteps N     total MD steps (default: 10000)\n"
        "  -vtk_dir DIR  VTK output directory (default: output)\n"
        "  -vtk_freq N   VTK output frequency (default: 1000)\n",
        prog);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Default configuration */
    CouplingState state;
    memset(&state, 0, sizeof(CouplingState));

    MPMConfig *cfg = &state.config;
    cfg->nelsx          = 20;
    cfg->nelsy          = 20;
    cfg->nelsz          = 20;
    cfg->sigma_bulk     = 1.0;
    cfg->sigma_defect   = 0.1;
    cfg->voltage_left   = 0.0;
    cfg->voltage_right  = 1.0;
    cfg->atom_charge    = 1.0;
    cfg->coupling_interval = 100;
    cfg->nr_max_iter    = 5;
    cfg->nr_tolerance   = 1e-9;
    cfg->mp_per_dir     = 1;   /* 1 MP per atom (not used in atom mode) */
    cfg->mp_type        = 2;   /* GIMP */

    int nsteps   = 10000;
    int vtk_freq = 1000;
    const char *input_script = NULL;
    const char *vtk_dir = "output";

    /* Parse command-line arguments */
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-nelsx") == 0 && i + 1 < argc) {
            cfg->nelsx = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nelsy") == 0 && i + 1 < argc) {
            cfg->nelsy = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nelsz") == 0 && i + 1 < argc) {
            cfg->nelsz = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-sigma_bulk") == 0 && i + 1 < argc) {
            cfg->sigma_bulk = atof(argv[++i]);
        } else if (strcmp(argv[i], "-sigma_def") == 0 && i + 1 < argc) {
            cfg->sigma_defect = atof(argv[++i]);
        } else if (strcmp(argv[i], "-vleft") == 0 && i + 1 < argc) {
            cfg->voltage_left = atof(argv[++i]);
        } else if (strcmp(argv[i], "-vright") == 0 && i + 1 < argc) {
            cfg->voltage_right = atof(argv[++i]);
        } else if (strcmp(argv[i], "-charge") == 0 && i + 1 < argc) {
            cfg->atom_charge = atof(argv[++i]);
        } else if (strcmp(argv[i], "-interval") == 0 && i + 1 < argc) {
            cfg->coupling_interval = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nsteps") == 0 && i + 1 < argc) {
            nsteps = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-vtk_dir") == 0 && i + 1 < argc) {
            vtk_dir = argv[++i];
        } else if (strcmp(argv[i], "-vtk_freq") == 0 && i + 1 < argc) {
            vtk_freq = atoi(argv[++i]);
        } else if (argv[i][0] != '-') {
            input_script = argv[i];
        } else {
            if (rank == 0) {
                fprintf(stderr, "Unknown option: %s\n", argv[i]);
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 1;
        }
        i++;
    }

    if (!input_script) {
        if (rank == 0) {
            fprintf(stderr, "Error: no LAMMPS input script specified\n");
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        printf("=== LAMMPS-MPM Electric Field Coupling ===\n");
        printf("Input script: %s\n", input_script);
        printf("Mesh: %d x %d x %d elements\n",
               cfg->nelsx, cfg->nelsy, cfg->nelsz);
        printf("Conductivity: bulk=%.4f, defect=%.4f\n",
               cfg->sigma_bulk, cfg->sigma_defect);
        printf("Voltage: left=%.4f, right=%.4f\n",
               cfg->voltage_left, cfg->voltage_right);
        printf("Charge: %.4f, Coupling interval: %d steps\n",
               cfg->atom_charge, cfg->coupling_interval);
        printf("Total steps: %d, VTK freq: %d\n", nsteps, vtk_freq);
    }

    /* Initialize LAMMPS and coupling */
    int ret = coupling_init(&state, argc, argv, input_script);
    if (ret != 0) {
        if (rank == 0)
            fprintf(stderr, "Coupling initialization failed\n");
        MPI_Finalize();
        return 1;
    }

    /* Write initial mesh VTK */
    if (rank == 0) {
        char fname[512];
        snprintf(fname, sizeof(fname), "%s/mesh.vtk", vtk_dir);
        vtk_write_mesh(&state.mesh, fname);
        printf("Wrote mesh VTK to %s\n", fname);
    }

    /* Run in chunks for VTK output */
    int steps_done = 0;
    int vtk_step = 0;
    while (steps_done < nsteps) {
        int chunk = vtk_freq;
        if (steps_done + chunk > nsteps)
            chunk = nsteps - steps_done;

        coupling_run(&state, chunk);
        steps_done += chunk;
        vtk_step++;

        /* Write VTK output */
        if (rank == 0 && state.mpdata.nmp > 0) {
            vtk_write_step(&state, vtk_dir, vtk_step);
            printf("Step %d/%d: wrote VTK (MPs=%d)\n",
                   steps_done, nsteps, state.mpdata.nmp);
        }
    }

    /* Cleanup */
    coupling_free(&state);

    if (rank == 0)
        printf("=== Coupling complete ===\n");

    MPI_Finalize();
    return 0;
}
