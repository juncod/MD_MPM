/*
 * mpm_electric.c — Newton-Raphson orchestrator for electric field solve
 * Ported from: ample_electric_3d.m (lines 49-76)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpm_electric.h"
#include "mpm_basis.h"
#include "mpm_assembly.h"
#include "mpm_solver.h"

/* --------------------------------------------------------------------------
 * Update MP potentials from converged nodal solution
 * Ported from updateMPs_Electric.m
 * --------------------------------------------------------------------------*/
void mpm_update_mp_potentials(MPData *mpdata, const double *phi_nodes)
{
    for (int i = 0; i < mpdata->nmp; i++) {
        MaterialPoint *mp = &mpdata->mp[i];
        double phi = 0.0;
        for (int j = 0; j < mp->nn; j++) {
            phi += mp->Svp[j] * phi_nodes[mp->nIN[j]];
        }
        mp->phi = phi;
    }
}

/* --------------------------------------------------------------------------
 * Full electric field solve
 *
 * Follows ample_electric_3d.m main loop:
 *   1. elemMPinfo → compute MP connectivity and basis functions
 *   2. detExtForce → assemble external current source
 *   3. Newton-Raphson loop:
 *      a. linSolve → solve for potential increment
 *      b. detMPs → assemble Ke and fint, compute gradPhi, J
 *      c. check convergence
 *   4. updateMPs → interpolate nodal phi to MPs
 * --------------------------------------------------------------------------*/
int mpm_electric_solve(CouplingState *state)
{
    Mesh        *mesh   = &state->mesh;
    MPData      *mpdata = &state->mpdata;
    SolverState *solver = &state->solver;
    MPMConfig   *cfg    = &state->config;

    int n = solver->n;
    int nD = 3;
    int NRitMax = cfg->nr_max_iter;
    double tol  = cfg->nr_tolerance;

    /* Step 1: Compute MP-element connectivity and basis functions */
    compute_mp_connectivity(mesh, mpdata, cfg);

    /* Step 2: Assemble external forces (Neumann BCs / volumetric sources) */
    assembly_compute_fext(solver, mpdata);

    /* Step 3: Identify free DOFs */
    assembly_compute_free_dofs(solver, mesh, mpdata);

    /* Step 4: Newton-Raphson iteration */
    double *phi_nodes = solver->phi;  /* solution vector */
    double *ddphi = (double *)mpm_calloc(n, sizeof(double));
    double *drct  = (double *)mpm_calloc(n, sizeof(double));

    /* Initialize out-of-balance force = fext (load factor = 1 for single lstp) */
    double *oobf = solver->rhs;
    memcpy(oobf, solver->fext, n * sizeof(double));

    memset(solver->frct, 0, n * sizeof(double));

    double fErr = 1.0;
    int NRit = 0;
    int converged = 0;

    while ((fErr > tol && NRit < NRitMax) || NRit < 2) {
        /* Solve for potential increment */
        /* Convert COO to CSR for CG solver */
        solver_coo_to_csr(solver);

        int ret = solver_solve(solver, mesh, ddphi, drct, NRit);
        if (ret != 0 && NRit > 0) {
            fprintf(stderr, "mpm_electric_solve: solver failed at NR iter %d\n", NRit);
            free(ddphi); free(drct);
            return -1;
        }

        /* Update nodal potential */
        for (int i = 0; i < n; i++) {
            phi_nodes[i] += ddphi[i];
            solver->frct[i] += drct[i];
        }

        /* Assemble conductance matrix and internal forces */
        solver_reset_coo(solver);
        assembly_compute_Ke_fint(solver, mpdata, phi_nodes, nD);

        /* Out-of-balance force: oobf = fext - fint + frct */
        double norm_oobf = 0.0;
        double norm_denom = 0.0;
        for (int i = 0; i < n; i++) {
            oobf[i] = solver->fext[i] - solver->fint[i] + solver->frct[i];
            norm_oobf += oobf[i] * oobf[i];
            double d = solver->fext[i] + solver->frct[i];
            norm_denom += d * d;
        }
        norm_oobf  = sqrt(norm_oobf);
        norm_denom = sqrt(norm_denom) + 1e-30; /* eps to avoid div-by-zero */

        fErr = norm_oobf / norm_denom;
        NRit++;

        fprintf(stdout, "  NR iteration %2d: error = %8.3e\n", NRit, fErr);

        if (fErr <= tol && NRit >= 2) {
            converged = 1;
            break;
        }
    }

    /* Step 5: Update MP potentials and electric fields */
    mpm_update_mp_potentials(mpdata, phi_nodes);

    /* Compute force on each MP: F = q * E */
    double q = cfg->atom_charge;
    for (int i = 0; i < mpdata->nmp; i++) {
        MaterialPoint *mp = &mpdata->mp[i];
        mp->force[0] = q * mp->E[0];
        mp->force[1] = q * mp->E[1];
        mp->force[2] = q * mp->E[2];
    }

    free(ddphi);
    free(drct);

    if (!converged) {
        fprintf(stderr, "mpm_electric_solve: NR did not converge (err=%e, tol=%e)\n",
                fErr, tol);
        return -1;
    }

    return 0;
}
