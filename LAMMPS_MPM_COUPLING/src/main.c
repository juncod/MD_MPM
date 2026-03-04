/*
 * main.c — LAMMPS-MPM Electric Field Coupling Driver
 *
 * All physics parameters are defined in the LAMMPS .in file as variables.
 * This driver simply:
 *   1. Opens LAMMPS
 *   2. Executes the .in file (setup + run 0)
 *   3. Reads variables, builds MPM mesh, registers callback
 *   4. Issues "run N" (N read from .in variable)
 *
 * Usage: mpirun -np N ./lammps_mpm_coupling input.in
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "mpm_types.h"
#include "lammps_callback.h"

/* LAMMPS library interface */
#include "library.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Find input script (first non-flag argument) */
    const char *input_script = NULL;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            input_script = argv[i];
            break;
        }
    }

    if (!input_script) {
        if (rank == 0) {
            fprintf(stderr,
                "Usage: mpirun -np N %s input.in\n\n"
                "All parameters are defined in the .in file as LAMMPS variables:\n"
                "  variable voltage_left   equal 0.0\n"
                "  variable voltage_right  equal 1.0\n"
                "  variable sigma_bulk     equal 1.0\n"
                "  variable sigma_defect   equal 0.1\n"
                "  variable atom_charge    equal 5.2e-6\n"
                "  variable target_h       equal 10.0\n"
                "  variable run_steps      equal 1000000\n"
                "  variable vtk_freq       equal 10000\n",
                argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        printf("=== LAMMPS-MPM Electric Field Coupling ===\n");
        printf("Input: %s\n", input_script);
    }

    /* Initialize coupling state */
    CouplingState state;
    memset(&state, 0, sizeof(CouplingState));

    /* Open LAMMPS */
    lammps_open(argc, argv, MPI_COMM_WORLD, &state.lmp);
    if (!state.lmp) {
        if (rank == 0) fprintf(stderr, "Failed to open LAMMPS\n");
        MPI_Finalize();
        return 1;
    }

    /* Execute .in, read variables, build mesh, register callback */
    int ret = coupling_init(&state, input_script);
    if (ret != 0) {
        if (rank == 0) fprintf(stderr, "Coupling initialization failed\n");
        coupling_free(&state);
        MPI_Finalize();
        return 1;
    }

    /* Run simulation (reads run_steps and vtk_freq from .in variables) */
    coupling_run(&state);

    /* Cleanup */
    coupling_free(&state);

    if (rank == 0)
        printf("=== Coupling complete ===\n");

    MPI_Finalize();
    return 0;
}
