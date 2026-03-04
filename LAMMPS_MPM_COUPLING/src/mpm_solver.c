/*
 * mpm_solver.c — Preconditioned Conjugate Gradient sparse solver
 * Ported from: linSolve_Electric.m
 *
 * Replaces UMFPACK with a pure-C PCG implementation.
 * K matrix is SPD for -div(σ∇φ)=0 → CG is guaranteed to converge.
 * Diagonal (Jacobi) preconditioner: M = diag(K).
 *
 * No external dependencies — only requires <math.h>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpm_solver.h"

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

    /* CG defaults */
    solver->cg_max_iter = 0;   /* 0 = auto (2*nfree) */
    solver->cg_tolerance = 1e-12;
}

/* --------------------------------------------------------------------------
 * COO → CSR conversion (sum duplicates)
 *
 * Algorithm:
 *   1. Count entries per row
 *   2. Build row pointers (prefix sum)
 *   3. Scatter entries into CSR arrays
 *   4. Sort each row by column index
 *   5. Merge duplicate (row,col) entries by summing values
 * --------------------------------------------------------------------------*/

/* Sort column indices within a row (insertion sort — rows are small) */
static void sort_row(int *cols, double *vals, int len)
{
    for (int i = 1; i < len; i++) {
        int    tc = cols[i];
        double tv = vals[i];
        int j = i - 1;
        while (j >= 0 && cols[j] > tc) {
            cols[j + 1] = cols[j];
            vals[j + 1] = vals[j];
            j--;
        }
        cols[j + 1] = tc;
        vals[j + 1] = tv;
    }
}

void solver_coo_to_csr(SolverState *solver)
{
    int n  = solver->n;
    int nz = solver->coo_nnz;

    /* Free old CSR */
    free(solver->Ap); solver->Ap = NULL;
    free(solver->Aj); solver->Aj = NULL;
    free(solver->Ax); solver->Ax = NULL;
    solver->nnz = 0;

    if (nz == 0) {
        solver->Ap = (int *)mpm_calloc(n + 1, sizeof(int));
        return;
    }

    /* Step 1: count entries per row */
    int *row_count = (int *)mpm_calloc(n, sizeof(int));
    for (int k = 0; k < nz; k++) {
        row_count[solver->coo_row[k]]++;
    }

    /* Step 2: build row pointers */
    int *Ap = (int *)mpm_calloc(n + 1, sizeof(int));
    for (int i = 0; i < n; i++) {
        Ap[i + 1] = Ap[i] + row_count[i];
    }

    /* Step 3: scatter COO into CSR (unsorted, with duplicates) */
    int *Aj = (int *)mpm_calloc(nz, sizeof(int));
    double *Ax = (double *)mpm_calloc(nz, sizeof(double));
    int *cursor = (int *)mpm_calloc(n, sizeof(int));
    memcpy(cursor, Ap, n * sizeof(int));

    for (int k = 0; k < nz; k++) {
        int r = solver->coo_row[k];
        int pos = cursor[r]++;
        Aj[pos] = solver->coo_col[k];
        Ax[pos] = solver->coo_val[k];
    }
    free(cursor);
    free(row_count);

    /* Step 4: sort each row by column index */
    for (int i = 0; i < n; i++) {
        int start = Ap[i], len = Ap[i + 1] - Ap[i];
        if (len > 1) {
            sort_row(&Aj[start], &Ax[start], len);
        }
    }

    /* Step 5: merge duplicates within each row */
    int new_nnz = 0;
    int *Ap_new = (int *)mpm_calloc(n + 1, sizeof(int));
    /* In-place compaction */
    for (int i = 0; i < n; i++) {
        int start = Ap[i], end = Ap[i + 1];
        Ap_new[i] = new_nnz;
        for (int k = start; k < end; ) {
            int col = Aj[k];
            double sum = 0.0;
            while (k < end && Aj[k] == col) {
                sum += Ax[k];
                k++;
            }
            Aj[new_nnz] = col;
            Ax[new_nnz] = sum;
            new_nnz++;
        }
    }
    Ap_new[n] = new_nnz;

    free(Ap);
    solver->Ap  = Ap_new;
    solver->Aj  = Aj;
    solver->Ax  = Ax;
    solver->nnz = new_nnz;
}

/* --------------------------------------------------------------------------
 * CSR matrix-vector multiply: y = A * x  (full matrix)
 * --------------------------------------------------------------------------*/
static void csr_matvec(const int *Ap, const int *Aj, const double *Ax,
                       int n, const double *x, double *y)
{
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = Ap[i]; k < Ap[i + 1]; k++) {
            s += Ax[k] * x[Aj[k]];
        }
        y[i] = s;
    }
}

/* --------------------------------------------------------------------------
 * Reduced CSR matvec: y = A(fd,fd) * x
 * perm[full_idx] = reduced_idx or -1
 * Only touches rows/cols in the free DOF set.
 * --------------------------------------------------------------------------*/
static void csr_matvec_reduced(const int *Ap, const int *Aj, const double *Ax,
                               int n_full, const int *perm,
                               int n_red, const double *x_red, double *y_red)
{
    /* Zero output */
    memset(y_red, 0, n_red * sizeof(double));

    for (int i = 0; i < n_full; i++) {
        int ri = perm[i];
        if (ri < 0) continue;
        double s = 0.0;
        for (int k = Ap[i]; k < Ap[i + 1]; k++) {
            int rj = perm[Aj[k]];
            if (rj >= 0) {
                s += Ax[k] * x_red[rj];
            }
        }
        y_red[ri] = s;
    }
}

/* --------------------------------------------------------------------------
 * Preconditioned Conjugate Gradient (Jacobi/diagonal preconditioner)
 *
 * Solves A_red * x = b for the free-DOF sub-system.
 * M^{-1} = 1/diag(A_red)
 *
 * Returns number of iterations (negative on failure).
 * --------------------------------------------------------------------------*/
static int pcg_solve_reduced(const int *Ap, const int *Aj, const double *Ax,
                             int n_full, const int *perm, const int *free_dofs,
                             int nfree, const double *b, double *x,
                             int max_iter, double tol)
{
    if (nfree == 0) return 0;

    /* Allocate CG work vectors */
    double *r  = (double *)mpm_calloc(nfree, sizeof(double));
    double *z  = (double *)mpm_calloc(nfree, sizeof(double));
    double *p  = (double *)mpm_calloc(nfree, sizeof(double));
    double *Ap_vec = (double *)mpm_calloc(nfree, sizeof(double));
    double *M_inv  = (double *)mpm_calloc(nfree, sizeof(double));

    /* Build diagonal preconditioner: M_inv[i] = 1 / A(fd[i], fd[i]) */
    for (int ii = 0; ii < nfree; ii++) {
        int i = free_dofs[ii];
        double diag = 0.0;
        for (int k = Ap[i]; k < Ap[i + 1]; k++) {
            if (Aj[k] == i) {
                diag = Ax[k];
                break;
            }
        }
        M_inv[ii] = (fabs(diag) > 1e-30) ? 1.0 / diag : 1.0;
    }

    /* r = b - A*x (initial residual) */
    csr_matvec_reduced(Ap, Aj, Ax, n_full, perm, nfree, x, Ap_vec);
    for (int i = 0; i < nfree; i++) {
        r[i] = b[i] - Ap_vec[i];
    }

    /* z = M^{-1} * r */
    double rz = 0.0;
    for (int i = 0; i < nfree; i++) {
        z[i] = M_inv[i] * r[i];
        rz += r[i] * z[i];
    }

    /* p = z */
    memcpy(p, z, nfree * sizeof(double));

    /* Compute initial residual norm for convergence check */
    double b_norm = 0.0;
    for (int i = 0; i < nfree; i++) b_norm += b[i] * b[i];
    b_norm = sqrt(b_norm);
    if (b_norm < 1e-30) b_norm = 1.0;

    double r_norm = 0.0;
    for (int i = 0; i < nfree; i++) r_norm += r[i] * r[i];
    r_norm = sqrt(r_norm);

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        if (r_norm / b_norm < tol) break;

        /* Ap_vec = A * p */
        csr_matvec_reduced(Ap, Aj, Ax, n_full, perm, nfree, p, Ap_vec);

        /* alpha = rz / (p' * Ap) */
        double pAp = 0.0;
        for (int i = 0; i < nfree; i++) pAp += p[i] * Ap_vec[i];

        if (fabs(pAp) < 1e-30) {
            /* Breakdown — matrix may be singular in this subspace */
            break;
        }
        double alpha = rz / pAp;

        /* x = x + alpha * p */
        /* r = r - alpha * Ap */
        for (int i = 0; i < nfree; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap_vec[i];
        }

        /* Recompute residual norm */
        r_norm = 0.0;
        for (int i = 0; i < nfree; i++) r_norm += r[i] * r[i];
        r_norm = sqrt(r_norm);

        /* z = M^{-1} * r */
        double rz_new = 0.0;
        for (int i = 0; i < nfree; i++) {
            z[i] = M_inv[i] * r[i];
            rz_new += r[i] * z[i];
        }

        /* beta = rz_new / rz */
        double beta = rz_new / (rz + 1e-30);
        rz = rz_new;

        /* p = z + beta * p */
        for (int i = 0; i < nfree; i++) {
            p[i] = z[i] + beta * p[i];
        }
    }

    free(r); free(z); free(p); free(Ap_vec); free(M_inv);

    return (r_norm / b_norm < tol) ? iter : -(iter + 1);
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
            ddphi[mesh->bc_node[i]] = mesh->bc_val[i];
            /* MATLAB: (1+sign(1-NRit))*vals = 1*vals at NRit=1 */
        }
    }

    int nfree = solver->nfree;
    if (nfree == 0) return 0;

    /* Build permutation: perm[full_idx] = reduced_idx or -1 */
    int *perm = solver->perm;
    memset(perm, -1, n * sizeof(int));
    for (int i = 0; i < nfree; i++) {
        perm[solver->free_dofs[i]] = i;
    }

    /* Compute K * ddphi_fixed (contribution from Dirichlet nodes) */
    double *Kd = (double *)mpm_calloc(n, sizeof(double));
    csr_matvec(solver->Ap, solver->Aj, solver->Ax, n, ddphi, Kd);

    /* Reduced RHS: b(fd) = oobf(fd) - Kd(fd) */
    double *b_red = (double *)mpm_calloc(nfree, sizeof(double));
    double *x_red = (double *)mpm_calloc(nfree, sizeof(double));

    for (int i = 0; i < nfree; i++) {
        b_red[i] = solver->rhs[solver->free_dofs[i]] - Kd[solver->free_dofs[i]];
        x_red[i] = 0.0;  /* initial guess = 0 */
    }

    /* Solve with PCG */
    int max_iter = solver->cg_max_iter > 0 ? solver->cg_max_iter : 2 * nfree;
    double tol = solver->cg_tolerance;

    int cg_iters = pcg_solve_reduced(
        solver->Ap, solver->Aj, solver->Ax,
        n, perm, solver->free_dofs, nfree,
        b_red, x_red, max_iter, tol);

    if (cg_iters < 0) {
        fprintf(stderr, "solver_solve: CG did not converge in %d iterations\n",
                -cg_iters);
        /* Continue with best approximation — don't fail hard */
    }

    /* Scatter solution back to full ddphi */
    for (int i = 0; i < nfree; i++) {
        ddphi[solver->free_dofs[i]] = x_red[i];
    }

    /* Compute reaction currents: drct(fixedNodes) = K * ddphi - oobf */
    memset(Kd, 0, n * sizeof(double));
    csr_matvec(solver->Ap, solver->Aj, solver->Ax, n, ddphi, Kd);

    for (int i = 0; i < mesh->nbc; i++) {
        int node = mesh->bc_node[i];
        drct[node] = Kd[node] - solver->rhs[node];
    }

    free(Kd);
    free(b_red);
    free(x_red);

    return 0;
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
    free(solver->Aj);        solver->Aj        = NULL;
    free(solver->Ax);        solver->Ax        = NULL;
}
