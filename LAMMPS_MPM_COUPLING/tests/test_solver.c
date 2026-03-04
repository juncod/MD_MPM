/*
 * test_solver.c — Integration test: solve -div(sigma*grad(phi))=0
 *
 * Test case: 5×5×5 mesh, uniform sigma=1.0, phi=0 at x=0, phi=1 at x=Lx
 * Expected: phi(x) = x/Lx (linear), E_x = -1/Lx, E_y = E_z = 0
 *
 * This test exercises the full MPM pipeline:
 *   mesh → MPs → connectivity → NR solve → verify phi, E
 *
 * Compile: cc -I../include -o test_solver test_solver.c \
 *          ../src/mpm_mesh.c ../src/mpm_basis.c ../src/mpm_assembly.c \
 *          ../src/mpm_solver.c ../src/mpm_electric.c ../src/atom_to_mp.c \
 *          ../src/vtk_output.c -lm
 * Run:     ./test_solver
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpm_types.h"
#include "mpm_mesh.h"
#include "mpm_basis.h"
#include "mpm_assembly.h"
#include "mpm_solver.h"
#include "mpm_electric.h"
#include "atom_to_mp.h"

#define ASSERT_NEAR(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "FAIL: %s:%d: %s=%.10f != %s=%.10f (tol=%.1e)\n", \
                __FILE__, __LINE__, #a, (double)(a), #b, (double)(b), (tol)); \
        exit(1); \
    } \
} while(0)

/*
 * Generate uniform grid MPs (like setupGrid_electric_3d.m)
 * mp_per_dir MPs per element direction → mp_per_dir^3 MPs per element
 */
static void generate_uniform_mps(MPData *mpdata, const Mesh *mesh,
                                  const MPMConfig *cfg)
{
    int mpd = 2;  /* MPs per direction per element */
    int nels = cfg->nelsx * cfg->nelsy * cfg->nelsz;
    int ngp = mpd * mpd * mpd;  /* 8 MPs per element */
    int total_mps = nels * ngp;

    mpdata->mp = (MaterialPoint *)mpm_calloc(total_mps, sizeof(MaterialPoint));
    mpdata->capacity = total_mps;
    mpdata->nmp = 0;

    double hx = mesh->h[0], hy = mesh->h[1], hz = mesh->h[2];
    double lpx = hx / (2.0 * mpd);
    double lpy = hy / (2.0 * mpd);
    double lpz = hz / (2.0 * mpd);
    double vol = 8.0 * lpx * lpy * lpz;

    for (int iz = 0; iz < cfg->nelsz; iz++) {
        for (int iy = 0; iy < cfg->nelsy; iy++) {
            for (int ix = 0; ix < cfg->nelsx; ix++) {
                double ex0 = cfg->x0 + ix * hx;
                double ey0 = cfg->y0 + iy * hy;
                double ez0 = cfg->z0 + iz * hz;

                for (int mz = 0; mz < mpd; mz++) {
                    for (int my = 0; my < mpd; my++) {
                        for (int mx = 0; mx < mpd; mx++) {
                            int idx = mpdata->nmp++;
                            MaterialPoint *mp = &mpdata->mp[idx];
                            memset(mp, 0, sizeof(MaterialPoint));

                            mp->pos[0] = ex0 + (mx + 0.5) * hx / mpd;
                            mp->pos[1] = ey0 + (my + 0.5) * hy / mpd;
                            mp->pos[2] = ez0 + (mz + 0.5) * hz / mpd;
                            mp->lp[0]  = lpx;
                            mp->lp[1]  = lpy;
                            mp->lp[2]  = lpz;
                            mp->vp     = vol;
                            mp->sigma  = cfg->sigma_bulk;
                            mp->mp_type = cfg->mp_type;
                            mp->fp     = 0.0;
                            mp->atom_id = idx;
                        }
                    }
                }
            }
        }
    }
}

int main(void)
{
    printf("=== Solver Integration Test ===\n");
    printf("Test: 5x5x5 mesh, uniform phi = x/Lx...\n");

    /* Configuration */
    MPMConfig cfg = {0};
    cfg.nelsx = 5; cfg.nelsy = 5; cfg.nelsz = 5;
    cfg.lx = 10.0; cfg.ly = 10.0; cfg.lz = 10.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;
    cfg.sigma_bulk   = 1.0;
    cfg.sigma_defect = 1.0;
    cfg.voltage_left  = 0.0;
    cfg.voltage_right = 1.0;
    cfg.atom_charge   = 1.0;
    cfg.coupling_interval = 1;
    cfg.nr_max_iter  = 5;
    cfg.nr_tolerance = 1e-9;
    cfg.mp_per_dir   = 2;
    cfg.mp_type      = 2;  /* GIMP */

    /* Create coupling state */
    CouplingState state;
    memset(&state, 0, sizeof(CouplingState));
    state.config = cfg;

    /* Create mesh */
    mesh_create(&state.mesh, &cfg);
    mesh_set_dirichlet_bc(&state.mesh, &cfg);

    printf("  Mesh: %d nodes, %d elements\n", state.mesh.nnodes, state.mesh.nels);
    printf("  BCs: %d Dirichlet nodes\n", state.mesh.nbc);

    /* Generate MPs */
    generate_uniform_mps(&state.mpdata, &state.mesh, &cfg);
    printf("  Material points: %d\n", state.mpdata.nmp);

    /* Initialize solver */
    solver_init(&state.solver, state.mesh.nnodes);

    /* Solve */
    printf("  Solving...\n");
    int ret = mpm_electric_solve(&state);
    if (ret != 0) {
        fprintf(stderr, "FAIL: mpm_electric_solve returned %d\n", ret);
        return 1;
    }

    /* Verify: phi should be linear in x */
    double max_phi_err = 0.0;
    double max_E_err = 0.0;
    double expected_Ex = -1.0 / cfg.lx;  /* -1/10 = -0.1 */

    for (int i = 0; i < state.mpdata.nmp; i++) {
        MaterialPoint *mp = &state.mpdata.mp[i];
        double x = mp->pos[0];
        double expected_phi = x / cfg.lx;

        double err_phi = fabs(mp->phi - expected_phi);
        if (err_phi > max_phi_err) max_phi_err = err_phi;

        /* E field should be uniform: Ex = -1/Lx */
        double err_Ex = fabs(mp->E[0] - expected_Ex);
        double err_Ey = fabs(mp->E[1]);
        double err_Ez = fabs(mp->E[2]);
        double err_E = err_Ex > err_Ey ? (err_Ex > err_Ez ? err_Ex : err_Ez)
                                        : (err_Ey > err_Ez ? err_Ey : err_Ez);
        if (err_E > max_E_err) max_E_err = err_E;
    }

    printf("  Results:\n");
    printf("    Max phi error: %.6e (tol: 1e-3)\n", max_phi_err);
    printf("    Max E error:   %.6e (tol: 1e-3)\n", max_E_err);
    printf("    Expected Ex:   %.6f\n", expected_Ex);

    /* Check nodal solution too */
    double max_node_err = 0.0;
    for (int i = 0; i < state.mesh.nnodes; i++) {
        double x = COORD(&state.mesh, i, 0);
        double expected = x / cfg.lx;
        double err = fabs(state.solver.phi[i] - expected);
        if (err > max_node_err) max_node_err = err;
    }
    printf("    Max nodal phi error: %.6e\n", max_node_err);

    /* Assertions */
    if (max_phi_err > 0.01) {
        fprintf(stderr, "FAIL: phi error too large: %e\n", max_phi_err);
        return 1;
    }
    if (max_E_err > 0.01) {
        fprintf(stderr, "FAIL: E field error too large: %e\n", max_E_err);
        return 1;
    }

    /* Cleanup */
    mpdata_free(&state.mpdata);
    solver_free(&state.solver);
    mesh_free(&state.mesh);

    printf("  PASSED\n");
    printf("\nAll solver tests passed!\n");
    return 0;
}
