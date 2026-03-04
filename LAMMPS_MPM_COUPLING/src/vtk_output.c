/*
 * vtk_output.c — VTK file output for debugging/visualization
 * Ported from: makeVtkMP_Electric_3D.m, postPro_Electric_3D.m
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vtk_output.h"

/* --------------------------------------------------------------------------
 * Write material point data to VTK file
 * --------------------------------------------------------------------------*/
void vtk_write_mp(const MPData *mpdata, const char *filename)
{
    FILE *fid = fopen(filename, "w");
    if (!fid) {
        fprintf(stderr, "vtk_write_mp: cannot open %s\n", filename);
        return;
    }

    int nmp = mpdata->nmp;

    fprintf(fid, "# vtk DataFile Version 2.0\n");
    fprintf(fid, "LAMMPS-MPM Electric field output\n");
    fprintf(fid, "ASCII\n");
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

    /* Points */
    fprintf(fid, "POINTS %d double\n", nmp);
    for (int i = 0; i < nmp; i++) {
        const MaterialPoint *mp = &mpdata->mp[i];
        fprintf(fid, "%f %f %f\n", mp->pos[0], mp->pos[1], mp->pos[2]);
    }
    fprintf(fid, "\n");

    /* Cells (each point is a vertex cell) */
    fprintf(fid, "CELLS %d %d\n", nmp, nmp * 2);
    for (int i = 0; i < nmp; i++) {
        fprintf(fid, "1 %d\n", i);
    }
    fprintf(fid, "\n");

    fprintf(fid, "CELL_TYPES %d\n", nmp);
    for (int i = 0; i < nmp; i++) {
        fprintf(fid, "1\n");  /* VTK_VERTEX */
    }
    fprintf(fid, "\n");

    /* Point data */
    fprintf(fid, "POINT_DATA %d\n", nmp);

    /* Electric potential (scalar) */
    fprintf(fid, "SCALARS phi FLOAT 1\n");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for (int i = 0; i < nmp; i++) {
        fprintf(fid, "%f\n", mpdata->mp[i].phi);
    }
    fprintf(fid, "\n");

    /* Conductivity (scalar) */
    fprintf(fid, "SCALARS sigma FLOAT 1\n");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for (int i = 0; i < nmp; i++) {
        fprintf(fid, "%f\n", mpdata->mp[i].sigma);
    }
    fprintf(fid, "\n");

    /* Electric field (vector) */
    fprintf(fid, "VECTORS E FLOAT\n");
    for (int i = 0; i < nmp; i++) {
        const MaterialPoint *mp = &mpdata->mp[i];
        fprintf(fid, "%f %f %f\n", mp->E[0], mp->E[1], mp->E[2]);
    }
    fprintf(fid, "\n");

    /* Force (vector) */
    fprintf(fid, "VECTORS force FLOAT\n");
    for (int i = 0; i < nmp; i++) {
        const MaterialPoint *mp = &mpdata->mp[i];
        fprintf(fid, "%f %f %f\n", mp->force[0], mp->force[1], mp->force[2]);
    }
    fprintf(fid, "\n");

    fclose(fid);
}

/* --------------------------------------------------------------------------
 * Write background mesh to VTK file (8-node hexahedra)
 * --------------------------------------------------------------------------*/
void vtk_write_mesh(const Mesh *mesh, const char *filename)
{
    FILE *fid = fopen(filename, "w");
    if (!fid) {
        fprintf(stderr, "vtk_write_mesh: cannot open %s\n", filename);
        return;
    }

    fprintf(fid, "# vtk DataFile Version 2.0\n");
    fprintf(fid, "MPM background mesh\n");
    fprintf(fid, "ASCII\n");
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

    /* Points */
    fprintf(fid, "POINTS %d double\n", mesh->nnodes);
    for (int i = 0; i < mesh->nnodes; i++) {
        fprintf(fid, "%f %f %f\n",
                mesh->coord[3*i], mesh->coord[3*i+1], mesh->coord[3*i+2]);
    }
    fprintf(fid, "\n");

    /* Cells (8-node hex) */
    fprintf(fid, "CELLS %d %d\n", mesh->nels, mesh->nels * 9);
    for (int e = 0; e < mesh->nels; e++) {
        fprintf(fid, "8");
        for (int n = 0; n < 8; n++) {
            fprintf(fid, " %d", mesh->etpl[8*e + n]);
        }
        fprintf(fid, "\n");
    }
    fprintf(fid, "\n");

    fprintf(fid, "CELL_TYPES %d\n", mesh->nels);
    for (int e = 0; e < mesh->nels; e++) {
        fprintf(fid, "12\n");  /* VTK_HEXAHEDRON */
    }
    fprintf(fid, "\n");

    fclose(fid);
}

/* --------------------------------------------------------------------------
 * Write step output (mesh + MP data)
 * --------------------------------------------------------------------------*/
void vtk_write_step(const CouplingState *state, const char *output_dir, int step)
{
    char fname[512];

    /* MP data */
    snprintf(fname, sizeof(fname), "%s/mpData_%d.vtk", output_dir, step);
    vtk_write_mp(&state->mpdata, fname);
}
