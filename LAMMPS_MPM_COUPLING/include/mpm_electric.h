/*
 * mpm_electric.h — Newton-Raphson orchestrator for electric field solve
 * Ported from: ample_electric_3d.m (NR loop)
 */
#ifndef MPM_ELECTRIC_H
#define MPM_ELECTRIC_H

#include "mpm_types.h"

/* Perform one complete electric field solve:
 * 1. Compute MP connectivity
 * 2. Assemble external forces
 * 3. Newton-Raphson iteration (assemble K, solve, update)
 * 4. Update MP potentials and compute E = -gradPhi
 *
 * Returns 0 on convergence, -1 on failure */
int mpm_electric_solve(CouplingState *state);

/* Update MP potentials from converged nodal solution */
void mpm_update_mp_potentials(MPData *mpdata, const double *phi_nodes);

#endif /* MPM_ELECTRIC_H */
