/*
 * test_mesh.c — Unit tests for 3D hex mesh generation
 *
 * Tests:
 *   1. 2×2×2 mesh: 27 nodes, 8 elements
 *   2. Node coordinates correctness
 *   3. Element topology (node ordering)
 *   4. Dirichlet BC assignment
 *
 * Compile: cc -I../include -o test_mesh test_mesh.c ../src/mpm_mesh.c -lm
 * Run:     ./test_mesh
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mpm_types.h"
#include "mpm_mesh.h"

#define ASSERT_EQ(a, b) do { \
    if ((a) != (b)) { \
        fprintf(stderr, "FAIL: %s:%d: %s=%d != %s=%d\n", \
                __FILE__, __LINE__, #a, (int)(a), #b, (int)(b)); \
        exit(1); \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "FAIL: %s:%d: %s=%.6f != %s=%.6f (tol=%.1e)\n", \
                __FILE__, __LINE__, #a, (double)(a), #b, (double)(b), (tol)); \
        exit(1); \
    } \
} while(0)

static void test_2x2x2_mesh(void)
{
    printf("Test: 2x2x2 mesh...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 2; cfg.nelsy = 2; cfg.nelsz = 2;
    cfg.lx = 4.0; cfg.ly = 6.0; cfg.lz = 8.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);

    /* Check counts */
    ASSERT_EQ(mesh.nnodes, 27);  /* (2+1)^3 = 27 */
    ASSERT_EQ(mesh.nels, 8);     /* 2^3 = 8 */

    /* Check element sizes */
    ASSERT_NEAR(mesh.h[0], 2.0, 1e-12);
    ASSERT_NEAR(mesh.h[1], 3.0, 1e-12);
    ASSERT_NEAR(mesh.h[2], 4.0, 1e-12);

    /* Check corner nodes */
    /* Node 0: (0,0,0) */
    ASSERT_NEAR(COORD(&mesh, 0, 0), 0.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 0, 1), 0.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 0, 2), 0.0, 1e-12);

    /* Last node: (4,6,8) */
    int last = mesh.nnodes - 1;
    ASSERT_NEAR(COORD(&mesh, last, 0), 4.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, last, 1), 6.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, last, 2), 8.0, 1e-12);

    /* Check first element (ix=0, iy=0, iz=0) */
    /* base=0, rowOff=3, planeOff=9 */
    ASSERT_EQ(ETPL(&mesh, 0, 0), 0);        /* N1 */
    ASSERT_EQ(ETPL(&mesh, 0, 1), 9);        /* N2: base+planeOff */
    ASSERT_EQ(ETPL(&mesh, 0, 2), 10);       /* N3: base+1+planeOff */
    ASSERT_EQ(ETPL(&mesh, 0, 3), 1);        /* N4: base+1 */
    ASSERT_EQ(ETPL(&mesh, 0, 4), 3);        /* N5: base+rowOff */
    ASSERT_EQ(ETPL(&mesh, 0, 5), 12);       /* N6: base+rowOff+planeOff */
    ASSERT_EQ(ETPL(&mesh, 0, 6), 13);       /* N7: base+rowOff+1+planeOff */
    ASSERT_EQ(ETPL(&mesh, 0, 7), 4);        /* N8: base+rowOff+1 */

    mesh_free(&mesh);
    printf("  PASSED\n");
}

static void test_dirichlet_bc(void)
{
    printf("Test: Dirichlet BCs...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 3; cfg.nelsy = 2; cfg.nelsz = 2;
    cfg.lx = 6.0; cfg.ly = 4.0; cfg.lz = 4.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;
    cfg.voltage_left  = 0.0;
    cfg.voltage_right = 5.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);
    mesh_set_dirichlet_bc(&mesh, &cfg);

    /* Left face nodes: x=0, (ny+1)*(nz+1) = 3*3 = 9 nodes */
    /* Right face nodes: x=6, 9 nodes */
    ASSERT_EQ(mesh.nbc, 18);  /* 9+9 */

    /* Check that all BC nodes have correct values */
    int n_left = 0, n_right = 0;
    for (int i = 0; i < mesh.nbc; i++) {
        double x = COORD(&mesh, mesh.bc_node[i], 0);
        if (fabs(x - 0.0) < 1e-10) {
            ASSERT_NEAR(mesh.bc_val[i], 0.0, 1e-12);
            n_left++;
        } else if (fabs(x - 6.0) < 1e-10) {
            ASSERT_NEAR(mesh.bc_val[i], 5.0, 1e-12);
            n_right++;
        } else {
            fprintf(stderr, "FAIL: BC node at unexpected x=%.2f\n", x);
            exit(1);
        }
    }
    ASSERT_EQ(n_left, 9);
    ASSERT_EQ(n_right, 9);

    mesh_free(&mesh);
    printf("  PASSED\n");
}

static void test_node_numbering(void)
{
    printf("Test: Node numbering order (x fastest)...\n");

    MPMConfig cfg = {0};
    cfg.nelsx = 3; cfg.nelsy = 2; cfg.nelsz = 2;
    cfg.lx = 3.0; cfg.ly = 2.0; cfg.lz = 2.0;
    cfg.x0 = 0.0; cfg.y0 = 0.0; cfg.z0 = 0.0;

    Mesh mesh = {0};
    mesh_create(&mesh, &cfg);

    /* Node numbering: x varies fastest
     * Node 0: (0,0,0), Node 1: (1,0,0), Node 2: (2,0,0), Node 3: (3,0,0)
     * Node 4: (0,1,0), ...  Node 7: (3,1,0)
     * Node 8: (0,2,0), ... Node 11: (3,2,0)
     * Node 12: (0,0,1), ... */
    ASSERT_NEAR(COORD(&mesh, 1, 0), 1.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 1, 1), 0.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 4, 0), 0.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 4, 1), 1.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 12, 0), 0.0, 1e-12);
    ASSERT_NEAR(COORD(&mesh, 12, 2), 1.0, 1e-12);

    mesh_free(&mesh);
    printf("  PASSED\n");
}

int main(void)
{
    printf("=== Mesh Unit Tests ===\n");
    test_2x2x2_mesh();
    test_dirichlet_bc();
    test_node_numbering();
    printf("\nAll mesh tests passed!\n");
    return 0;
}
