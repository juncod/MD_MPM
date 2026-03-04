/*
 * atom_to_mp.h — Atom ↔ Material Point conversion
 */
#ifndef ATOM_TO_MP_H
#define ATOM_TO_MP_H

#include "mpm_types.h"

/* Convert atom positions + CNA data to material points
 * positions: natoms×3 array (row-major)
 * cna:       natoms array of CNA values (5=FCC bulk, other=defect)
 * ids:       natoms array of atom IDs
 */
void atoms_to_mps(MPData *mpdata, const MPMConfig *cfg,
                  const double *positions, const int *cna,
                  const int *ids, int natoms);

/* Extract forces from converged MP solution back to atom array
 * forces_out: natoms×3 (row-major), indexed by atom order
 */
void mps_to_forces(const MPData *mpdata, double *forces_out, int natoms);

/* Free MP connectivity data (nIN, Svp, dSvp arrays) */
void mpdata_free_connectivity(MPData *mpdata);

/* Free all MP data */
void mpdata_free(MPData *mpdata);

#endif /* ATOM_TO_MP_H */
