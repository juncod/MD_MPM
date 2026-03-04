/*
 * mpm_basis.h — GIMP/MPM shape functions and MP-element connectivity
 * Ported from: SvpGIMP.m, SvpMPM.m, MPMbasis.m, elemForMP.m, nodesForMP.m
 */
#ifndef MPM_BASIS_H
#define MPM_BASIS_H

#include "mpm_types.h"

/* 1D GIMP basis function (5 piecewise regions) */
void basis_gimp_1d(double xp, double xv, double h, double lp,
                   double *Svp, double *dSvp);

/* 1D standard MPM basis function (piecewise linear) */
void basis_mpm_1d(double xp, double xv, double h,
                  double *Svp, double *dSvp);

/* 3D basis function via tensor product of 1D functions
 * Returns Svp (scalar) and dSvp[3] (gradient) */
void basis_3d(const Mesh *mesh, const MaterialPoint *mp, int node,
              double *Svp, double dSvp[3]);

/* Find element index containing position pos[3] using O(1) grid lookup
 * Returns element index (0-based), or -1 if outside domain */
int elem_for_pos(const Mesh *mesh, const MPMConfig *cfg, const double pos[3]);

/* Find all elements overlapping MP's GIMP domain
 * Returns number of elements found; elem_list must have enough capacity */
int elems_for_mp(const Mesh *mesh, const MPMConfig *cfg,
                 const double pos[3], const double lp[3],
                 int *elem_list, int max_elems);

/* Extract unique sorted nodes from a set of elements
 * Returns number of unique nodes; node_list must have enough capacity */
int nodes_for_elems(const Mesh *mesh, const int *elem_list, int nelem,
                    int *node_list, int max_nodes);

/* Compute connectivity and basis functions for all MPs */
void compute_mp_connectivity(Mesh *mesh, MPData *mpdata, const MPMConfig *cfg);

#endif /* MPM_BASIS_H */
