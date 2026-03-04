/*
 * test_basis.c — Unit tests for GIMP/MPM shape functions
 *
 * Tests:
 *   1. GIMP 1D: 5-region piecewise values
 *   2. Partition of unity: sum(N_i) = 1 at any point
 *   3. Linear field gradient: if phi = x, then gradPhi = [1,0,0]
 *   4. Element lookup correctness
 *
 * Compile: cc -I../include -o test_basis test_basis.c ../src/mpm_basis.c ../src/mpm_mesh.c -lm
 * Run:     ./test_basis
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpm_types.h"
#include "mpm_mesh.h"
#include "mpm_basis.h"

#define ASSERT_NEAR(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "FAIL: %s:%d: %s=%.10f != %s=%.10f (tol=%.1e)\n", \
                __FILE__, __LINE__, #a, (double)(a), #b, (double)(b), (tol)); \
        exit(1); \
    } \
} while(0)

static void test_gimp_1d_regions(void)
{
    printf("Test: GIMP 1D regions...\n");

    double h = 1.0, lp = 0.25;

    /* Test at node position (xp = xv): region C, should be ~1 */
    double S, dS;
    basis_gimp_1d(5.0, 5.0, h, lp, &S, &dS);
    double expected = 1.0 - (0.0 + lp*lp) / (2.0*h*lp);
    ASSERT_NEAR(S, expected, 1e-12);   /* 1 - lp/(2h) = 0.875 */
    ASSERT_NEAR(dS, 0.0, 1e-12);       /* derivative = 0 at center */

    /* Test outside range: should be 0 */
    basis_gimp_1d(5.0, 8.0, h, lp, &S, &dS);
    ASSERT_NEAR(S, 0.0, 1e-12);
    ASSERT_NEAR(dS, 0.0, 1e-12);

    /* Test region A: d = xp - xv in (-h-lp, -h+lp] */
    /* xp - xv = -h (midpoint of region A) */
    basis_gimp_1d(5.0, 6.0, h, lp, &S, &dS);
    /* d = -1.0, region B: 1 + d/h = 0 */
    /* Actually d = -1.0 is boundary of A and B. Check region B edge: */
    basis_gimp_1d(5.0, 5.75, h, lp, &S, &dS);
    /* d = -0.75, region B: d in (-h+lp, -lp] = (-0.75, -0.25] */
    ASSERT_NEAR(S, 1.0 + (-0.75)/h, 1e-12);  /* 0.25 */
    ASSERT_NEAR(dS, 1.0/h, 1e-12);

    printf("  PASSED\n");
}

static void test_partition_of_unity(void)
{
    printf("Test: Partition of unity (sum N_i = 1)...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 4; cfg.nelsy = 4; cfg.nelsz = 4;
    cfg.lx = 4.0; cfg.ly = 4.0; cfg.lz = 4.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);

    /* Test at several interior points */
    double test_points[][3] = {
        {2.0, 2.0, 2.0},   /* center */
        {1.3, 2.7, 0.8},   /* arbitrary */
        {0.5, 0.5, 0.5},   /* near corner */
        {3.5, 3.5, 3.5},   /* near opposite corner */
    };
    int npts = 4;

    for (int p = 0; p < npts; p++) {
        MaterialPoint mp = {0};
        mp.pos[0] = test_points[p][0];
        mp.pos[1] = test_points[p][1];
        mp.pos[2] = test_points[p][2];
        mp.lp[0] = mesh.h[0] / 4.0;  /* GIMP half-length */
        mp.lp[1] = mesh.h[1] / 4.0;
        mp.lp[2] = mesh.h[2] / 4.0;
        mp.mp_type = 2;  /* GIMP */

        /* Find overlapping elements */
        int elem_buf[64];
        int nelem = elems_for_mp(&mesh, &cfg, mp.pos, mp.lp, elem_buf, 64);

        /* Find unique nodes */
        int node_buf[256];
        int nn = nodes_for_elems(&mesh, elem_buf, nelem, node_buf, 256);

        /* Sum basis functions */
        double sum_S = 0.0;
        double sum_dS[3] = {0.0, 0.0, 0.0};
        for (int j = 0; j < nn; j++) {
            double S, dS[3];
            basis_3d(&mesh, &mp, node_buf[j], &S, dS);
            sum_S += S;
            sum_dS[0] += dS[0];
            sum_dS[1] += dS[1];
            sum_dS[2] += dS[2];
        }

        ASSERT_NEAR(sum_S, 1.0, 1e-10);
        ASSERT_NEAR(sum_dS[0], 0.0, 1e-10);
        ASSERT_NEAR(sum_dS[1], 0.0, 1e-10);
        ASSERT_NEAR(sum_dS[2], 0.0, 1e-10);
    }

    mesh_free(&mesh);
    printf("  PASSED\n");
}

static void test_linear_field_gradient(void)
{
    printf("Test: Linear field gradient (phi=x → gradPhi=[1,0,0])...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 5; cfg.nelsy = 5; cfg.nelsz = 5;
    cfg.lx = 10.0; cfg.ly = 10.0; cfg.lz = 10.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);

    /* Set phi_nodes = x coordinate at each node */
    double *phi = (double *)mpm_calloc(mesh.nnodes, sizeof(double));
    for (int i = 0; i < mesh.nnodes; i++) {
        phi[i] = COORD(&mesh, i, 0);  /* phi = x */
    }

    /* Test gradient at interior point */
    MaterialPoint mp = {0};
    mp.pos[0] = 4.3; mp.pos[1] = 6.1; mp.pos[2] = 3.7;
    mp.lp[0] = mesh.h[0] / 4.0;
    mp.lp[1] = mesh.h[1] / 4.0;
    mp.lp[2] = mesh.h[2] / 4.0;
    mp.mp_type = 2;

    int elem_buf[64];
    int nelem = elems_for_mp(&mesh, &cfg, mp.pos, mp.lp, elem_buf, 64);
    int node_buf[256];
    int nn = nodes_for_elems(&mesh, elem_buf, nelem, node_buf, 256);

    /* Compute gradient: gradPhi = sum_j dS_j * phi_j */
    double gradPhi[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < nn; j++) {
        double S, dS[3];
        basis_3d(&mesh, &mp, node_buf[j], &S, dS);
        gradPhi[0] += dS[0] * phi[node_buf[j]];
        gradPhi[1] += dS[1] * phi[node_buf[j]];
        gradPhi[2] += dS[2] * phi[node_buf[j]];
    }

    ASSERT_NEAR(gradPhi[0], 1.0, 1e-10);
    ASSERT_NEAR(gradPhi[1], 0.0, 1e-10);
    ASSERT_NEAR(gradPhi[2], 0.0, 1e-10);

    free(phi);
    mesh_free(&mesh);
    printf("  PASSED\n");
}

static void test_elem_lookup(void)
{
    printf("Test: Element lookup...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 5; cfg.nelsy = 4; cfg.nelsz = 3;
    cfg.lx = 10.0; cfg.ly = 8.0; cfg.lz = 6.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);

    /* Point in element (2,1,1): ix=2, iy=1, iz=1
     * Expected element index: 1*5*4 + 1*5 + 2 = 27 */
    double pos[3] = {5.5, 3.5, 3.5};  /* center of element (2,1,1) */
    int elem = elem_for_pos(&mesh, &cfg, pos);

    /* Element index: iz*nx*ny + iy*nx + ix = 1*20 + 1*5 + 2 = 27 */
    if (elem != 27) {
        fprintf(stderr, "FAIL: elem_for_pos returned %d, expected 27\n", elem);
        exit(1);
    }

    /* Point outside domain */
    double pos_out[3] = {-1.0, 0.0, 0.0};
    int elem_out = elem_for_pos(&mesh, &cfg, pos_out);
    if (elem_out != -1) {
        fprintf(stderr, "FAIL: expected -1 for out-of-domain point, got %d\n", elem_out);
        exit(1);
    }

    mesh_free(&mesh);
    printf("  PASSED\n");
}

int main(void)
{
    printf("=== Basis Function Unit Tests ===\n");
    test_gimp_1d_regions();
    test_partition_of_unity();
    test_linear_field_gradient();
    test_elem_lookup();
    printf("\nAll basis tests passed!\n");
    return 0;
}
