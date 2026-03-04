/*
 * mpm_assembly.h — Stiffness matrix and force vector assembly
 * Ported from: detMPs_Electric.m, detExtForce_Electric.m, detFDoFs_Electric.m
 */
#ifndef MPM_ASSEMBLY_H
#define MPM_ASSEMBLY_H

#include "mpm_types.h"

/* Assemble global conductance matrix (COO) and internal force vector
 * Also computes gradPhi, J, E at each MP */
void assembly_compute_Ke_fint(SolverState *solver, MPData *mpdata,
                              const double *phi_nodes, int nD);

/* Assemble external current source vector from MP volumetric sources */
void assembly_compute_fext(SolverState *solver, const MPData *mpdata);

/* Identify free DOFs (active nodes minus Dirichlet BC nodes) */
void assembly_compute_free_dofs(SolverState *solver, const Mesh *mesh,
                                const MPData *mpdata);

#endif /* MPM_ASSEMBLY_H */
