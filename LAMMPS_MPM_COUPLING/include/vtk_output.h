/*
 * vtk_output.h — VTK file output for debugging/visualization
 * Ported from: makeVtkMP_Electric_3D.m, postPro_Electric_3D.m
 */
#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include "mpm_types.h"

/* Write material point data to VTK file
 * Outputs: phi (scalar), J (vector), sigma (scalar) */
void vtk_write_mp(const MPData *mpdata, const char *filename);

/* Write background mesh to VTK file (hex cells) */
void vtk_write_mesh(const Mesh *mesh, const char *filename);

/* Write both mesh and MP data for a given timestep
 * Files written to: output_dir/mesh.vtk, output_dir/mpData_<step>.vtk */
void vtk_write_step(const CouplingState *state, const char *output_dir, int step);

#endif /* VTK_OUTPUT_H */
