/*
 * lammps_callback.h — LAMMPS fix external callback for MPM coupling
 */
#ifndef LAMMPS_CALLBACK_H
#define LAMMPS_CALLBACK_H

#include "mpm_types.h"

/* Callback function signature for fix external pf/callback
 * Called by LAMMPS at each timestep where fix external is active */
void mpm_force_callback(void *caller, int64_t timestep,
                        int nlocal, int *ids, double **x, double **f);

/* Initialize LAMMPS coupling: open LAMMPS, register callback */
int coupling_init(CouplingState *state, int argc, char **argv,
                  const char *input_script);

/* Run LAMMPS for nsteps */
int coupling_run(CouplingState *state, int nsteps);

/* Cleanup */
void coupling_free(CouplingState *state);

#endif /* LAMMPS_CALLBACK_H */
