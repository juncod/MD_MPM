/*
 * mpm_basis.c — GIMP/MPM shape functions and MP-element connectivity
 * Ported from: SvpGIMP.m, SvpMPM.m, MPMbasis.m, elemForMP.m, nodesForMP.m,
 *              elemMPinfo_Electric.m
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpm_basis.h"

/* --------------------------------------------------------------------------
 * 1D GIMP basis function — 5 piecewise regions (A–E)
 * Ported from SvpGIMP.m
 * --------------------------------------------------------------------------*/
void basis_gimp_1d(double xp, double xv, double h, double lp,
                   double *Svp, double *dSvp)
{
    double d = xp - xv;  /* signed distance */

    if (d > (-h - lp) && d <= (-h + lp)) {
        /* A: partial overlap, left boundary */
        double t = h + lp + d;
        *Svp  = (t * t) / (4.0 * h * lp);
        *dSvp = t / (2.0 * h * lp);
    }
    else if (d > (-h + lp) && d <= (-lp)) {
        /* B: full overlap, left element */
        *Svp  = 1.0 + d / h;
        *dSvp = 1.0 / h;
    }
    else if (d > (-lp) && d <= lp) {
        /* C: partial overlap, straddling two elements */
        *Svp  = 1.0 - (d * d + lp * lp) / (2.0 * h * lp);
        *dSvp = -d / (h * lp);
    }
    else if (d > lp && d <= (h - lp)) {
        /* D: full overlap, right element */
        *Svp  = 1.0 - d / h;
        *dSvp = -1.0 / h;
    }
    else if (d > (h - lp) && d <= (h + lp)) {
        /* E: partial overlap, right boundary */
        double t = h + lp - d;
        *Svp  = (t * t) / (4.0 * h * lp);
        *dSvp = -t / (2.0 * h * lp);
    }
    else {
        /* No overlap */
        *Svp  = 0.0;
        *dSvp = 0.0;
    }
}

/* --------------------------------------------------------------------------
 * 1D standard MPM basis function (piecewise linear)
 * Ported from SvpMPM.m
 * --------------------------------------------------------------------------*/
void basis_mpm_1d(double xp, double xv, double h,
                  double *Svp, double *dSvp)
{
    double d = xp - xv;

    if (d > -h && d <= 0.0) {
        *Svp  = 1.0 + d / h;
        *dSvp = 1.0 / h;
    }
    else if (d > 0.0 && d <= h) {
        *Svp  = 1.0 - d / h;
        *dSvp = -1.0 / h;
    }
    else {
        *Svp  = 0.0;
        *dSvp = 0.0;
    }
}

/* --------------------------------------------------------------------------
 * 3D basis function via tensor product
 * Ported from MPMbasis.m
 * --------------------------------------------------------------------------*/
void basis_3d(const Mesh *mesh, const MaterialPoint *mp, int node,
              double *Svp_out, double dSvp_out[3])
{
    double S[3], dS[3];

    for (int i = 0; i < 3; i++) {
        double xp = mp->pos[i];
        double xv = mesh->coord[3 * node + i];
        double h  = mesh->h[i];

        if (mp->mp_type == 2) {
            basis_gimp_1d(xp, xv, h, mp->lp[i], &S[i], &dS[i]);
        } else {
            basis_mpm_1d(xp, xv, h, &S[i], &dS[i]);
        }
    }

    /* Tensor product: Svp = S[0]*S[1]*S[2] */
    *Svp_out = S[0] * S[1] * S[2];

    /* Gradient: dSvp[i] = dS[i] * prod(S[j], j!=i)
     * 3D index pattern from MPMbasis.m: indx = [2 3; 1 3; 1 2] */
    dSvp_out[0] = dS[0] * S[1] * S[2];
    dSvp_out[1] = dS[1] * S[0] * S[2];
    dSvp_out[2] = dS[2] * S[0] * S[1];
}

/* --------------------------------------------------------------------------
 * O(1) element lookup for a position
 * --------------------------------------------------------------------------*/
int elem_for_pos(const Mesh *mesh, const MPMConfig *cfg, const double pos[3])
{
    int nx = cfg->nelsx, ny = cfg->nelsy, nz = cfg->nelsz;
    double x0 = cfg->x0, y0 = cfg->y0, z0 = cfg->z0;

    int ix = (int)((pos[0] - x0) / mesh->h[0]);
    int iy = (int)((pos[1] - y0) / mesh->h[1]);
    int iz = (int)((pos[2] - z0) / mesh->h[2]);

    /* Clamp to valid range */
    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
        return -1;

    return iz * (nx * ny) + iy * nx + ix;
}

/* --------------------------------------------------------------------------
 * Find all elements overlapping MP's GIMP domain (AABB test)
 * Ported from elemForMP.m but using grid indexing for O(1) per dimension
 * --------------------------------------------------------------------------*/
int elems_for_mp(const Mesh *mesh, const MPMConfig *cfg,
                 const double pos[3], const double lp[3],
                 int *elem_list, int max_elems)
{
    int nx = cfg->nelsx, ny = cfg->nelsy, nz = cfg->nelsz;
    double x0 = cfg->x0, y0 = cfg->y0, z0 = cfg->z0;

    /* Compute element index ranges that overlap [pos-lp, pos+lp] */
    int ix_lo = (int)((pos[0] - lp[0] - x0) / mesh->h[0]);
    int ix_hi = (int)((pos[0] + lp[0] - x0) / mesh->h[0]);
    int iy_lo = (int)((pos[1] - lp[1] - y0) / mesh->h[1]);
    int iy_hi = (int)((pos[1] + lp[1] - y0) / mesh->h[1]);
    int iz_lo = (int)((pos[2] - lp[2] - z0) / mesh->h[2]);
    int iz_hi = (int)((pos[2] + lp[2] - z0) / mesh->h[2]);

    /* Clamp */
    if (ix_lo < 0)  ix_lo = 0;  if (ix_hi >= nx) ix_hi = nx - 1;
    if (iy_lo < 0)  iy_lo = 0;  if (iy_hi >= ny) iy_hi = ny - 1;
    if (iz_lo < 0)  iz_lo = 0;  if (iz_hi >= nz) iz_hi = nz - 1;

    int count = 0;
    for (int iz = iz_lo; iz <= iz_hi && count < max_elems; iz++) {
        for (int iy = iy_lo; iy <= iy_hi && count < max_elems; iy++) {
            for (int ix = ix_lo; ix <= ix_hi && count < max_elems; ix++) {
                elem_list[count++] = iz * (nx * ny) + iy * nx + ix;
            }
        }
    }
    return count;
}

/* --------------------------------------------------------------------------
 * Extract unique sorted nodes from element list
 * Ported from nodesForMP.m
 * --------------------------------------------------------------------------*/
static int int_compare(const void *a, const void *b)
{
    return (*(const int *)a) - (*(const int *)b);
}

int nodes_for_elems(const Mesh *mesh, const int *elem_list, int nelem,
                    int *node_list, int max_nodes)
{
    /* Collect all nodes (with duplicates) */
    int total = nelem * 8;
    int *buf = (int *)mpm_calloc(total, sizeof(int));
    int count = 0;
    for (int e = 0; e < nelem; e++) {
        for (int n = 0; n < 8; n++) {
            buf[count++] = mesh->etpl[8 * elem_list[e] + n];
        }
    }

    /* Sort and deduplicate */
    qsort(buf, count, sizeof(int), int_compare);

    int unique = 0;
    for (int i = 0; i < count && unique < max_nodes; i++) {
        if (i == 0 || buf[i] != buf[i - 1]) {
            node_list[unique++] = buf[i];
        }
    }

    free(buf);
    return unique;
}

/* --------------------------------------------------------------------------
 * Compute connectivity and basis functions for all MPs
 * Ported from elemMPinfo_Electric.m
 * --------------------------------------------------------------------------*/
void compute_mp_connectivity(Mesh *mesh, MPData *mpdata, const MPMConfig *cfg)
{
    /* Reset active elements */
    memset(mesh->eInA, 0, mesh->nels * sizeof(int));

    /* Temporary buffers (max possible: 3^3=27 elements, 4^3=64 nodes for GIMP) */
    int elem_buf[64];
    int node_buf[256];

    for (int i = 0; i < mpdata->nmp; i++) {
        MaterialPoint *mp = &mpdata->mp[i];

        /* Find overlapping elements */
        int nelem = elems_for_mp(mesh, cfg, mp->pos, mp->lp,
                                 elem_buf, 64);

        /* Mark active elements */
        for (int e = 0; e < nelem; e++) {
            mesh->eInA[elem_buf[e]] = 1;
        }

        /* Find unique nodes */
        int nn = nodes_for_elems(mesh, elem_buf, nelem, node_buf, 256);

        /* Free old connectivity if any */
        free(mp->nIN);
        free(mp->Svp);
        free(mp->dSvp);

        /* Allocate new connectivity */
        mp->nn   = nn;
        mp->nIN  = (int *)mpm_calloc(nn, sizeof(int));
        mp->Svp  = (double *)mpm_calloc(nn, sizeof(double));
        mp->dSvp = (double *)mpm_calloc(3 * nn, sizeof(double));
        mp->nSMe = nn * nn;

        memcpy(mp->nIN, node_buf, nn * sizeof(int));

        /* Compute basis functions at each connected node */
        for (int j = 0; j < nn; j++) {
            double S, dS[3];
            basis_3d(mesh, mp, mp->nIN[j], &S, dS);
            mp->Svp[j] = S;
            mp->dSvp[3 * j + 0] = dS[0];
            mp->dSvp[3 * j + 1] = dS[1];
            mp->dSvp[3 * j + 2] = dS[2];
        }
    }
}
