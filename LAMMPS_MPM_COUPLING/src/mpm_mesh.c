/*
 * mpm_mesh.c — 3D hexahedral mesh generation
 * Ported from: MPM_EXAMPLE_ELECTRIC_3D/setup/formCoord3D.m
 *
 * Node numbering: x fastest, then y, then z
 * Element topology: 8-node hex with node ordering matching formCoord3D.m
 *   N1(x-,y-,z-), N2(x-,y-,z+), N3(x+,y-,z+), N4(x+,y-,z-),
 *   N5(x-,y+,z-), N6(x-,y+,z+), N7(x+,y+,z+), N8(x+,y+,z-)
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpm_mesh.h"

void mesh_create(Mesh *mesh, const MPMConfig *cfg)
{
    int nx = cfg->nelsx, ny = cfg->nelsy, nz = cfg->nelsz;
    double lx = cfg->lx, ly = cfg->ly, lz = cfg->lz;
    double x0 = cfg->x0, y0 = cfg->y0, z0 = cfg->z0;

    mesh->nnodes = (nx + 1) * (ny + 1) * (nz + 1);
    mesh->nels   = nx * ny * nz;
    mesh->h[0]   = lx / nx;
    mesh->h[1]   = ly / ny;
    mesh->h[2]   = lz / nz;

    /* Allocate arrays */
    mesh->coord = (double *)mpm_calloc(mesh->nnodes * 3, sizeof(double));
    mesh->etpl  = (int *)mpm_calloc(mesh->nels * 8, sizeof(int));
    mesh->eInA  = (int *)mpm_calloc(mesh->nels, sizeof(int));

    /* Node generation: loop z → y → x (x fastest) */
    int node = 0;
    for (int k = 0; k <= nz; k++) {
        double z = z0 + lz * k / nz;
        for (int j = 0; j <= ny; j++) {
            double y = y0 + ly * j / ny;
            for (int i = 0; i <= nx; i++) {
                double x = x0 + lx * i / nx;
                mesh->coord[3 * node + 0] = x;
                mesh->coord[3 * node + 1] = y;
                mesh->coord[3 * node + 2] = z;
                node++;
            }
        }
    }

    /* Element generation (0-based node indices) */
    int rowOff   = nx + 1;
    int planeOff = (nx + 1) * (ny + 1);
    int nel = 0;
    for (int iz = 0; iz < nz; iz++) {
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                int base = iz * planeOff + iy * rowOff + ix;
                int *e = &mesh->etpl[8 * nel];
                e[0] = base;                             /* N1: (x-,y-,z-) */
                e[1] = base + planeOff;                  /* N2: (x-,y-,z+) */
                e[2] = base + 1 + planeOff;              /* N3: (x+,y-,z+) */
                e[3] = base + 1;                         /* N4: (x+,y-,z-) */
                e[4] = base + rowOff;                    /* N5: (x-,y+,z-) */
                e[5] = base + rowOff + planeOff;         /* N6: (x-,y+,z+) */
                e[6] = base + rowOff + 1 + planeOff;     /* N7: (x+,y+,z+) */
                e[7] = base + rowOff + 1;                /* N8: (x+,y+,z-) */
                nel++;
            }
        }
    }

    /* BC arrays initially empty */
    mesh->nbc     = 0;
    mesh->bc_node = NULL;
    mesh->bc_val  = NULL;
}

void mesh_set_dirichlet_bc(Mesh *mesh, const MPMConfig *cfg)
{
    /* Count BC nodes: left face (x=x0) and right face (x=x0+lx) */
    double x_left  = cfg->x0;
    double x_right = cfg->x0 + cfg->lx;
    double tol     = mesh->h[0] * 1e-6;

    int count = 0;
    for (int i = 0; i < mesh->nnodes; i++) {
        double x = mesh->coord[3 * i];
        if (fabs(x - x_left) < tol || fabs(x - x_right) < tol)
            count++;
    }

    /* Free old BCs if any */
    free(mesh->bc_node);
    free(mesh->bc_val);

    mesh->nbc     = count;
    mesh->bc_node = (int *)mpm_calloc(count, sizeof(int));
    mesh->bc_val  = (double *)mpm_calloc(count, sizeof(double));

    int idx = 0;
    for (int i = 0; i < mesh->nnodes; i++) {
        double x = mesh->coord[3 * i];
        if (fabs(x - x_left) < tol) {
            mesh->bc_node[idx] = i;
            mesh->bc_val[idx]  = cfg->voltage_left;
            idx++;
        } else if (fabs(x - x_right) < tol) {
            mesh->bc_node[idx] = i;
            mesh->bc_val[idx]  = cfg->voltage_right;
            idx++;
        }
    }
}

void mesh_free(Mesh *mesh)
{
    free(mesh->coord);   mesh->coord   = NULL;
    free(mesh->etpl);    mesh->etpl    = NULL;
    free(mesh->eInA);    mesh->eInA    = NULL;
    free(mesh->bc_node); mesh->bc_node = NULL;
    free(mesh->bc_val);  mesh->bc_val  = NULL;
    mesh->nnodes = 0;
    mesh->nels   = 0;
    mesh->nbc    = 0;
}
