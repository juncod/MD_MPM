/*
 * mpm_assembly.c — Stiffness matrix and force vector assembly
 * Ported from: detMPs_Electric.m, detExtForce_Electric.m, detFDoFs_Electric.m
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpm_assembly.h"

/* --------------------------------------------------------------------------
 * Assemble global conductance matrix (COO) and internal force vector
 * Ported from detMPs_Electric.m
 *
 * For each MP:
 *   G = dNx (3 × nn)
 *   D = sigma * I (3 × 3)
 *   Ke += G' * D * G * vol     (nn × nn contribution)
 *   fint += G' * D * G * phi_nodes * vol
 *   gradPhi = G * phi_nodes
 *   J = -D * gradPhi
 * --------------------------------------------------------------------------*/
void assembly_compute_Ke_fint(SolverState *solver, MPData *mpdata,
                              const double *phi_nodes, int nD)
{
    /* Reset internal force */
    memset(solver->fint, 0, solver->n * sizeof(double));

    /* Count total COO entries needed */
    int total_entries = 0;
    for (int i = 0; i < mpdata->nmp; i++) {
        total_entries += mpdata->mp[i].nSMe;
    }

    /* Ensure COO buffer capacity */
    if (total_entries > solver->coo_capacity) {
        solver->coo_capacity = total_entries + total_entries / 4;
        solver->coo_row = (int *)mpm_realloc(solver->coo_row,
                          solver->coo_capacity * sizeof(int));
        solver->coo_col = (int *)mpm_realloc(solver->coo_col,
                          solver->coo_capacity * sizeof(int));
        solver->coo_val = (double *)mpm_realloc(solver->coo_val,
                          solver->coo_capacity * sizeof(double));
    }
    solver->coo_nnz = 0;

    for (int mi = 0; mi < mpdata->nmp; mi++) {
        MaterialPoint *mp = &mpdata->mp[mi];
        int nn = mp->nn;
        if (nn == 0) continue;

        double sigma  = mp->sigma;
        double vol    = mp->vp;

        /* G = dSvp (3 × nn), stored row-major in mp->dSvp */
        /* Compute gradPhi = G * phi_nodes(nIN) — (3×1) */
        double gradPhi[3] = {0.0, 0.0, 0.0};
        for (int j = 0; j < nn; j++) {
            double phi_j = phi_nodes[mp->nIN[j]];
            gradPhi[0] += mp->dSvp[3 * j + 0] * phi_j;
            gradPhi[1] += mp->dSvp[3 * j + 1] * phi_j;
            gradPhi[2] += mp->dSvp[3 * j + 2] * phi_j;
        }

        /* D = sigma * I, so D*gradPhi = sigma * gradPhi */
        double flux[3];
        flux[0] = sigma * gradPhi[0];
        flux[1] = sigma * gradPhi[1];
        flux[2] = sigma * gradPhi[2];

        /* Store MP fields */
        mp->gradPhi[0] = gradPhi[0];
        mp->gradPhi[1] = gradPhi[1];
        mp->gradPhi[2] = gradPhi[2];
        mp->E[0] = -gradPhi[0];        /* E = -gradPhi */
        mp->E[1] = -gradPhi[1];
        mp->E[2] = -gradPhi[2];

        /* Ke contribution: kp(a,b) = sum_d (G(d,a) * sigma * G(d,b)) * vol
         * fint contribution: fp(a) = sum_d (G(d,a) * flux(d)) * vol */
        int base = solver->coo_nnz;
        for (int a = 0; a < nn; a++) {
            /* fint(nIN[a]) += G'(a,:) * flux * vol */
            double fp_a = 0.0;
            for (int d = 0; d < nD; d++) {
                fp_a += mp->dSvp[3 * a + d] * flux[d];
            }
            solver->fint[mp->nIN[a]] += fp_a * vol;

            /* Ke entries for row a */
            for (int b = 0; b < nn; b++) {
                double kp_ab = 0.0;
                for (int d = 0; d < nD; d++) {
                    kp_ab += mp->dSvp[3 * a + d] * sigma * mp->dSvp[3 * b + d];
                }
                solver->coo_row[base] = mp->nIN[a];
                solver->coo_col[base] = mp->nIN[b];
                solver->coo_val[base] = kp_ab * vol;
                base++;
            }
        }
        solver->coo_nnz = base;
    }
}

/* --------------------------------------------------------------------------
 * Assemble external current source vector
 * Ported from detExtForce_Electric.m
 * --------------------------------------------------------------------------*/
void assembly_compute_fext(SolverState *solver, const MPData *mpdata)
{
    memset(solver->fext, 0, solver->n * sizeof(double));

    for (int i = 0; i < mpdata->nmp; i++) {
        const MaterialPoint *mp = &mpdata->mp[i];
        double fp = mp->fp;
        if (fp == 0.0) continue;

        for (int j = 0; j < mp->nn; j++) {
            solver->fext[mp->nIN[j]] += fp * mp->Svp[j];
        }
    }
}

/* --------------------------------------------------------------------------
 * Identify free DOFs
 * Ported from detFDoFs_Electric.m
 * --------------------------------------------------------------------------*/
void assembly_compute_free_dofs(SolverState *solver, const Mesh *mesh,
                                const MPData *mpdata)
{
    int n = solver->n;

    /* Mark active nodes (in active elements) */
    char *active = (char *)mpm_calloc(n, sizeof(char));
    for (int e = 0; e < mesh->nels; e++) {
        if (mesh->eInA[e]) {
            for (int k = 0; k < 8; k++) {
                active[mesh->etpl[8 * e + k]] = 1;
            }
        }
    }

    /* Mark Dirichlet BC nodes as not free */
    char *is_bc = (char *)mpm_calloc(n, sizeof(char));
    for (int i = 0; i < mesh->nbc; i++) {
        is_bc[mesh->bc_node[i]] = 1;
    }

    /* Collect free DOFs: active AND not BC */
    free(solver->free_dofs);
    solver->free_dofs = (int *)mpm_calloc(n, sizeof(int));
    solver->nfree = 0;

    for (int i = 0; i < n; i++) {
        if (active[i] && !is_bc[i]) {
            solver->free_dofs[solver->nfree++] = i;
        }
    }

    free(active);
    free(is_bc);
}
