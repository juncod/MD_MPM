/*
 * lammps_callback.h — LAMMPS fix external callback for MPM coupling
 *
 * All physics parameters are read from LAMMPS variables defined in .in file.
 * Callback frequency is controlled by fix external Nevery parameter.
 */
#ifndef LAMMPS_CALLBACK_H
#define LAMMPS_CALLBACK_H

#include "mpm_types.h"

/* Callback function for fix external pf/callback
 * Called by LAMMPS every Nevery steps (set in .in file) */
void mpm_force_callback(void *caller, int64_t timestep,
                        int nlocal, int *ids, double **x, double **f);

/* Initialize: execute .in, read variables, build mesh, register callback */
int coupling_init(CouplingState *state, const char *input_script);

/* Run: read run_steps from .in variable, issue "run N" */
int coupling_run(CouplingState *state);

/* Cleanup */
void coupling_free(CouplingState *state);

#endif /* LAMMPS_CALLBACK_H */
