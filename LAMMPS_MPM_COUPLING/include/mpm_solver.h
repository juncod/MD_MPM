/*
 * mpm_solver.h — UMFPACK sparse direct solver wrapper
 * Ported from: linSolve_Electric.m
 */
#ifndef MPM_SOLVER_H
#define MPM_SOLVER_H

#include "mpm_types.h"

/* Initialize solver state for n DOFs */
void solver_init(SolverState *solver, int n);

/* Convert COO triplets to CSC format (summing duplicates) */
void solver_coo_to_csc(SolverState *solver);

/* Perform UMFPACK symbolic + numeric factorization */
int solver_factorize(SolverState *solver);

/* Solve the reduced system for free DOFs with Dirichlet BCs
 * ddphi: potential increment (output, size n)
 * drct:  reaction currents (output, size n)
 * NRit:  Newton-Raphson iteration number (1-based)
 * Returns 0 on success */
int solver_solve(SolverState *solver, const Mesh *mesh,
                 double *ddphi, double *drct, int NRit);

/* Free solver memory */
void solver_free(SolverState *solver);

/* Reset COO assembly (clear triplets, keep allocated memory) */
void solver_reset_coo(SolverState *solver);

#endif /* MPM_SOLVER_H */
