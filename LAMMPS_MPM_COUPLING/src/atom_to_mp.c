/*
 * atom_to_mp.c — Atom ↔ Material Point conversion
 *
 * Converts LAMMPS atom positions + CNA data to MPM material points.
 * CNA=5 (FCC bulk) → sigma_bulk, others → sigma_defect
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom_to_mp.h"

/* --------------------------------------------------------------------------
 * Convert atoms to material points
 *
 * Each atom becomes one MP. Position is mapped from LAMMPS box coords
 * to the MPM mesh coordinate system. Volume is uniformly distributed.
 * Conductivity is assigned based on CNA value.
 * --------------------------------------------------------------------------*/
void atoms_to_mps(MPData *mpdata, const MPMConfig *cfg,
                  const double *positions, const int *cna,
                  const int *ids, int natoms)
{
    /* Ensure capacity */
    if (mpdata->capacity < natoms) {
        mpdata->mp = (MaterialPoint *)mpm_realloc(
            mpdata->mp, natoms * sizeof(MaterialPoint));
        mpdata->capacity = natoms;
    }
    mpdata->nmp = natoms;

    /* Compute uniform volume per atom */
    double total_volume = cfg->lx * cfg->ly * cfg->lz;
    double vol_per_atom = total_volume / natoms;

    /* GIMP half-lengths based on Cu lattice constant (~3.615 Å) */
    double a_Cu = 3.615;
    double lp_default = a_Cu / 2.0;  /* ~1.8 Å */

    /* Clamp half-lengths to half element size */
    double lp[3];
    for (int d = 0; d < 3; d++) {
        lp[d] = lp_default;
        double max_lp = cfg->mp_type == 2 ? 0.49 * (d == 0 ? cfg->lx / cfg->nelsx :
                        d == 1 ? cfg->ly / cfg->nelsy : cfg->lz / cfg->nelsz) : 0.0;
        if (lp[d] > max_lp) lp[d] = max_lp;
    }

    for (int i = 0; i < natoms; i++) {
        MaterialPoint *mp = &mpdata->mp[i];
        memset(mp, 0, sizeof(MaterialPoint));

        /* Position: LAMMPS coords are already in the simulation box
         * Map to MPM mesh: pos = atom_pos (assuming box origin = mesh origin) */
        mp->pos[0] = positions[3 * i + 0];
        mp->pos[1] = positions[3 * i + 1];
        mp->pos[2] = positions[3 * i + 2];

        /* GIMP domain half-lengths */
        mp->lp[0] = lp[0];
        mp->lp[1] = lp[1];
        mp->lp[2] = lp[2];

        /* Volume */
        mp->vp = vol_per_atom;

        /* Conductivity based on CNA */
        if (cna && cna[i] == 5) {
            mp->sigma = cfg->sigma_bulk;    /* FCC bulk */
        } else {
            mp->sigma = cfg->sigma_defect;  /* Defect site */
        }

        /* No external current source */
        mp->fp = 0.0;

        /* Atom ID */
        mp->atom_id = ids ? ids[i] : i;

        /* MP type */
        mp->mp_type = cfg->mp_type;

        /* Connectivity will be set by compute_mp_connectivity */
        mp->nn   = 0;
        mp->nIN  = NULL;
        mp->Svp  = NULL;
        mp->dSvp = NULL;
    }
}

/* --------------------------------------------------------------------------
 * Extract forces from MP solution back to atom force array
 * --------------------------------------------------------------------------*/
void mps_to_forces(const MPData *mpdata, double *forces_out, int natoms)
{
    /* Zero output */
    memset(forces_out, 0, 3 * natoms * sizeof(double));

    /* Each MP maps 1:1 to an atom (by index order) */
    int n = mpdata->nmp < natoms ? mpdata->nmp : natoms;
    for (int i = 0; i < n; i++) {
        const MaterialPoint *mp = &mpdata->mp[i];
        forces_out[3 * i + 0] = mp->force[0];
        forces_out[3 * i + 1] = mp->force[1];
        forces_out[3 * i + 2] = mp->force[2];
    }
}

/* --------------------------------------------------------------------------
 * Free MP connectivity data
 * --------------------------------------------------------------------------*/
void mpdata_free_connectivity(MPData *mpdata)
{
    for (int i = 0; i < mpdata->nmp; i++) {
        free(mpdata->mp[i].nIN);  mpdata->mp[i].nIN  = NULL;
        free(mpdata->mp[i].Svp);  mpdata->mp[i].Svp  = NULL;
        free(mpdata->mp[i].dSvp); mpdata->mp[i].dSvp = NULL;
        mpdata->mp[i].nn = 0;
    }
}

/* --------------------------------------------------------------------------
 * Free all MP data
 * --------------------------------------------------------------------------*/
void mpdata_free(MPData *mpdata)
{
    if (mpdata->mp) {
        mpdata_free_connectivity(mpdata);
        free(mpdata->mp);
        mpdata->mp = NULL;
    }
    mpdata->nmp = 0;
    mpdata->capacity = 0;
}
