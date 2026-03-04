/*
 * mpm_mesh.h — 3D hexahedral mesh generation
 * Ported from: MPM_EXAMPLE_ELECTRIC_3D/setup/formCoord3D.m
 */
#ifndef MPM_MESH_H
#define MPM_MESH_H

#include "mpm_types.h"

/* Create 3D hex mesh with nelsx×nelsy×nelsz elements over [x0,x0+lx]×[y0,y0+ly]×[z0,z0+lz] */
void mesh_create(Mesh *mesh, const MPMConfig *cfg);

/* Set Dirichlet BCs: phi=voltage_left at x=x0, phi=voltage_right at x=x0+lx */
void mesh_set_dirichlet_bc(Mesh *mesh, const MPMConfig *cfg);

/* Free all mesh memory */
void mesh_free(Mesh *mesh);

#endif /* MPM_MESH_H */
