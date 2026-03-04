/*
 * mpm_solver.c — UMFPACK sparse direct solver wrapper
 * Ported from: linSolve_Electric.m
 *
 * Uses SuiteSparse UMFPACK for sparse LU factorization.
 * Dirichlet BCs applied by extracting the free-DOF sub-system.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpm_solver.h"
#include <umfpack.h>

/* --------------------------------------------------------------------------
 * Initialize solver state
 * --------------------------------------------------------------------------*/
void solver_init(SolverState *solver, int n)
{
    memset(solver, 0, sizeof(SolverState));
    solver->n = n;

    solver->phi  = (double *)mpm_calloc(n, sizeof(double));
    solver->rhs  = (double *)mpm_calloc(n, sizeof(double));
    solver->fint = (double *)mpm_calloc(n, sizeof(double));
    solver->fext = (double *)mpm_calloc(n, sizeof(double));
    solver->frct = (double *)mpm_calloc(n, sizeof(double));
    solver->perm = (int *)mpm_calloc(n, sizeof(int));

    /* Initial COO capacity */
    solver->coo_capacity = 1024;
    solver->coo_nnz = 0;
    solver->coo_row = (int *)mpm_calloc(solver->coo_capacity, sizeof(int));
    solver->coo_col = (int *)mpm_calloc(solver->coo_capacity, sizeof(int));
    solver->coo_val = (double *)mpm_calloc(solver->coo_capacity, sizeof(double));
}

/* --------------------------------------------------------------------------
 * Convert COO to CSC format using UMFPACK triplet_to_col
 * --------------------------------------------------------------------------*/
void solver_coo_to_csc(SolverState *solver)
{
    int n = solver->n;
    int nz = solver->coo_nnz;

    /* Free old CSC */
    free(solver->Ap);
    free(solver->Ai);
    free(solver->Ax);

    if (nz == 0) {
        solver->nnz = 0;
        solver->Ap = (int *)mpm_calloc(n + 1, sizeof(int));
        solver->Ai = NULL;
        solver->Ax = NULL;
        return;
    }

    /* Allocate CSC arrays (max size = nz, may be smaller after summing) */
    solver->Ap = (int *)mpm_calloc(n + 1, sizeof(int));
    solver->Ai = (int *)mpm_calloc(nz, sizeof(int));
    solver->Ax = (double *)mpm_calloc(nz, sizeof(double));

    int status = umfpack_di_triplet_to_col(
        n, n, nz,
        solver->coo_row, solver->coo_col, solver->coo_val,
        solver->Ap, solver->Ai, solver->Ax, NULL);

    if (status != UMFPACK_OK) {
        fprintf(stderr, "umfpack_di_triplet_to_col failed: %d\n", status);
        return;
    }

    solver->nnz = solver->Ap[n];
}

/* --------------------------------------------------------------------------
 * UMFPACK factorize
 * --------------------------------------------------------------------------*/
int solver_factorize(SolverState *solver)
{
    int n = solver->n;

    /* Free old factorizations */
    if (solver->Symbolic) {
        umfpack_di_free_symbolic(&solver->Symbolic);
        solver->Symbolic = NULL;
    }
    if (solver->Numeric) {
        umfpack_di_free_numeric(&solver->Numeric);
        solver->Numeric = NULL;
    }

    int status = umfpack_di_symbolic(n, n,
        solver->Ap, solver->Ai, solver->Ax,
        &solver->Symbolic, NULL, NULL);
    if (status != UMFPACK_OK) {
        fprintf(stderr, "umfpack_di_symbolic failed: %d\n", status);
        return -1;
    }

    status = umfpack_di_numeric(
        solver->Ap, solver->Ai, solver->Ax,
        solver->Symbolic, &solver->Numeric, NULL, NULL);
    if (status != UMFPACK_OK) {
        fprintf(stderr, "umfpack_di_numeric failed: %d\n", status);
        return -1;
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Solve for potential increment and reaction currents
 * Ported from linSolve_Electric.m
 *
 * NRit: Newton-Raphson iteration (1-based)
 *   NRit=0: skip solve (need K first)
 *   NRit=1: apply Dirichlet BCs to ddphi, then solve
 *   NRit>1: only solve for free DOFs (BCs already applied)
 * --------------------------------------------------------------------------*/
int solver_solve(SolverState *solver, const Mesh *mesh,
                 double *ddphi, double *drct, int NRit)
{
    int n = solver->n;
    memset(ddphi, 0, n * sizeof(double));
    memset(drct, 0, n * sizeof(double));

    if (NRit == 0) return 0;

    /* Apply Dirichlet BCs at NRit=1 */
    if (NRit == 1) {
        for (int i = 0; i < mesh->nbc; i++) {
            ddphi[mesh->bc_node[i]] = 2.0 * mesh->bc_val[i];
            /* Factor of 2 matches MATLAB: (1+sign(1-NRit))*fixedValues = 2*vals at NRit=1 */
        }
    }
    /* At NRit>1, ddphi at fixed nodes stays 0 (already converged) */

    int nfree = solver->nfree;
    if (nfree == 0) return 0;

    /* Build reduced RHS: b_free = oobf(fd) - Ke(fd, fixedNodes) * ddphi(fixedNodes)
     * oobf is stored in solver->rhs
     * We need to compute K * ddphi for the fixed nodes contribution */

    /* First, compute K * ddphi (full vector) — only nonzero for fixed nodes */
    double *Kd = (double *)mpm_calloc(n, sizeof(double));
    /* CSC multiply: for each column j, add Ax[k]*ddphi[j] to Kd[Ai[k]] */
    for (int j = 0; j < n; j++) {
        if (ddphi[j] == 0.0) continue;
        for (int k = solver->Ap[j]; k < solver->Ap[j + 1]; k++) {
            Kd[solver->Ai[k]] += solver->Ax[k] * ddphi[j];
        }
    }

    /* Reduced RHS: b(fd) = oobf(fd) - Kd(fd) */
    double *b_red = (double *)mpm_calloc(nfree, sizeof(double));
    for (int i = 0; i < nfree; i++) {
        b_red[i] = solver->rhs[solver->free_dofs[i]] - Kd[solver->free_dofs[i]];
    }

    /* Extract reduced K matrix (free DOFs only) in CSC format */
    /* Build permutation: perm[full_idx] = reduced_idx or -1 */
    int *perm = solver->perm;
    memset(perm, -1, n * sizeof(int));
    for (int i = 0; i < nfree; i++) {
        perm[solver->free_dofs[i]] = i;
    }

    /* Count non-zeros in reduced matrix */
    int nnz_red = 0;
    for (int j = 0; j < n; j++) {
        if (perm[j] < 0) continue;
        for (int k = solver->Ap[j]; k < solver->Ap[j + 1]; k++) {
            if (perm[solver->Ai[k]] >= 0) nnz_red++;
        }
    }

    int *Ap_red = (int *)mpm_calloc(nfree + 1, sizeof(int));
    int *Ai_red = (int *)mpm_calloc(nnz_red, sizeof(int));
    double *Ax_red = (double *)mpm_calloc(nnz_red, sizeof(double));

    /* Fill reduced CSC */
    int ptr = 0;
    for (int jj = 0; jj < nfree; jj++) {
        int j = solver->free_dofs[jj];
        Ap_red[jj] = ptr;
        for (int k = solver->Ap[j]; k < solver->Ap[j + 1]; k++) {
            int ri = perm[solver->Ai[k]];
            if (ri >= 0) {
                Ai_red[ptr] = ri;
                Ax_red[ptr] = solver->Ax[k];
                ptr++;
            }
        }
    }
    Ap_red[nfree] = ptr;

    /* Solve reduced system with UMFPACK */
    void *Sym_red = NULL, *Num_red = NULL;
    double *x_red = (double *)mpm_calloc(nfree, sizeof(double));
    int status;

    status = umfpack_di_symbolic(nfree, nfree, Ap_red, Ai_red, Ax_red,
                                 &Sym_red, NULL, NULL);
    if (status != UMFPACK_OK) {
        fprintf(stderr, "solver_solve: symbolic failed (%d)\n", status);
        goto cleanup;
    }

    status = umfpack_di_numeric(Ap_red, Ai_red, Ax_red,
                                Sym_red, &Num_red, NULL, NULL);
    if (status != UMFPACK_OK) {
        fprintf(stderr, "solver_solve: numeric failed (%d)\n", status);
        goto cleanup;
    }

    status = umfpack_di_solve(UMFPACK_A, Ap_red, Ai_red, Ax_red,
                              x_red, b_red, Num_red, NULL, NULL);
    if (status != UMFPACK_OK) {
        fprintf(stderr, "solver_solve: solve failed (%d)\n", status);
        goto cleanup;
    }

    /* Scatter solution back to full ddphi */
    for (int i = 0; i < nfree; i++) {
        ddphi[solver->free_dofs[i]] = x_red[i];
    }

    /* Compute reaction currents: drct(fixedNodes) = K * ddphi - oobf */
    /* Recompute K * ddphi with the full solution */
    memset(Kd, 0, n * sizeof(double));
    for (int j = 0; j < n; j++) {
        if (ddphi[j] == 0.0) continue;
        for (int k = solver->Ap[j]; k < solver->Ap[j + 1]; k++) {
            Kd[solver->Ai[k]] += solver->Ax[k] * ddphi[j];
        }
    }
    for (int i = 0; i < mesh->nbc; i++) {
        int node = mesh->bc_node[i];
        drct[node] = Kd[node] - solver->rhs[node];
    }

cleanup:
    umfpack_di_free_symbolic(&Sym_red);
    umfpack_di_free_numeric(&Num_red);
    free(Ap_red); free(Ai_red); free(Ax_red);
    free(x_red); free(b_red); free(Kd);

    return (status == UMFPACK_OK) ? 0 : -1;
}

/* --------------------------------------------------------------------------
 * Reset COO assembly
 * --------------------------------------------------------------------------*/
void solver_reset_coo(SolverState *solver)
{
    solver->coo_nnz = 0;
}

/* --------------------------------------------------------------------------
 * Free solver
 * --------------------------------------------------------------------------*/
void solver_free(SolverState *solver)
{
    free(solver->phi);       solver->phi       = NULL;
    free(solver->rhs);       solver->rhs       = NULL;
    free(solver->fint);      solver->fint      = NULL;
    free(solver->fext);      solver->fext      = NULL;
    free(solver->frct);      solver->frct      = NULL;
    free(solver->perm);      solver->perm      = NULL;
    free(solver->free_dofs); solver->free_dofs = NULL;
    free(solver->coo_row);   solver->coo_row   = NULL;
    free(solver->coo_col);   solver->coo_col   = NULL;
    free(solver->coo_val);   solver->coo_val   = NULL;
    free(solver->Ap);        solver->Ap        = NULL;
    free(solver->Ai);        solver->Ai        = NULL;
    free(solver->Ax);        solver->Ax        = NULL;
    if (solver->Symbolic) umfpack_di_free_symbolic(&solver->Symbolic);
    if (solver->Numeric)  umfpack_di_free_numeric(&solver->Numeric);
}
