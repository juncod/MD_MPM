/*
 * mpm_types.h — Core data structures for LAMMPS-MPM Electric Field Coupling
 *
 * All struct definitions used across the MPM solver modules.
 */
#ifndef MPM_TYPES_H
#define MPM_TYPES_H

#include <stdlib.h>

/* ---------------------------------------------------------------------------
 * Configuration
 * -------------------------------------------------------------------------*/
typedef struct {
    double lx, ly, lz;              /* domain dimensions (Angstrom) */
    int    nelsx, nelsy, nelsz;     /* element counts per direction */
    double sigma_bulk;              /* conductivity for bulk (CNA=5) */
    double sigma_defect;            /* conductivity for defect sites */
    double voltage_left;            /* Dirichlet BC at x=x0           */
    double voltage_right;           /* Dirichlet BC at x=x0+lx        */
    double atom_charge;             /* F = q*E charge parameter        */
    int    nr_max_iter;             /* Newton-Raphson max iterations    */
    double nr_tolerance;            /* NR convergence tolerance         */
    double x0, y0, z0;             /* domain origin                    */
    double target_h;                /* target element edge length (Å)   */
    int    mp_per_dir;              /* material points per element dir  */
    int    mp_type;                 /* 1=MPM, 2=GIMP                   */
} MPMConfig;

/* ---------------------------------------------------------------------------
 * 3D Hexahedral Mesh
 * -------------------------------------------------------------------------*/
typedef struct {
    int    nnodes;                  /* total node count */
    int    nels;                    /* total element count */
    double h[3];                    /* element size per direction */
    double *coord;                  /* nnodes × 3 (row-major) */
    int    *etpl;                   /* nels × 8   (row-major, 0-based) */
    int    nbc;                     /* number of Dirichlet BC nodes */
    int    *bc_node;                /* bc_node[nbc] — node indices */
    double *bc_val;                 /* bc_val[nbc]  — prescribed values */
    int    *eInA;                   /* nels × 1 — active element flag */
} Mesh;

/* ---------------------------------------------------------------------------
 * Material Point
 * -------------------------------------------------------------------------*/
typedef struct {
    double pos[3];                  /* coordinates */
    double lp[3];                   /* GIMP domain half-lengths */
    double vp;                      /* current volume */
    double sigma;                   /* conductivity */
    double phi;                     /* electric potential */
    double gradPhi[3];              /* potential gradient */
    double E[3];                    /* electric field = -gradPhi */
    double force[3];                /* electromigration force = q*E */
    double fp;                      /* volumetric current source */
    /* connectivity (variable length, allocated per MP) */
    int    nn;                      /* number of connected nodes */
    int    *nIN;                    /* node indices [nn] */
    double *Svp;                    /* basis function values [nn] */
    double *dSvp;                   /* basis gradients [3*nn] row-major */
    int    nSMe;                    /* nn^2 — stiffness entries */
    int    atom_id;                 /* owning atom ID (-1 if none) */
    int    mp_type;                 /* 1=MPM, 2=GIMP */
} MaterialPoint;

/* ---------------------------------------------------------------------------
 * Material Point Data (collection)
 * -------------------------------------------------------------------------*/
typedef struct {
    int            nmp;             /* number of material points */
    int            capacity;        /* allocated capacity */
    MaterialPoint *mp;              /* array of material points */
} MPData;

/* ---------------------------------------------------------------------------
 * Sparse Solver State (COO → CSR + Conjugate Gradient)
 * -------------------------------------------------------------------------*/
typedef struct {
    int    n;                       /* matrix dimension (nnodes) */
    /* COO triplets (for assembly) */
    int    coo_capacity;
    int    coo_nnz;
    int    *coo_row;
    int    *coo_col;
    double *coo_val;
    /* CSR format (for CG solver) */
    int    nnz;
    int    *Ap;                     /* row pointers [n+1] */
    int    *Aj;                     /* column indices [nnz] */
    double *Ax;                     /* values [nnz] */
    /* CG solver parameters */
    int    cg_max_iter;             /* max CG iterations (default: 2*n) */
    double cg_tolerance;            /* CG convergence tolerance */
    /* Vectors */
    double *phi;                    /* solution vector [n] */
    double *rhs;                    /* right-hand side [n] */
    double *fint;                   /* internal force [n] */
    double *fext;                   /* external force [n] */
    double *frct;                   /* reaction currents [n] */
    /* Free DOFs */
    int    *free_dofs;
    int    nfree;
    /* Dirichlet-reduced system */
    int    *perm;                   /* full→reduced index map [n] */
} SolverState;

/* ---------------------------------------------------------------------------
 * Coupling State (top-level, owns everything)
 * -------------------------------------------------------------------------*/
typedef struct {
    MPMConfig   config;
    Mesh        mesh;
    MPData      mpdata;
    SolverState solver;
    double      *prev_forces;       /* reserved for future caching */
    int         prev_forces_size;
    void        *lmp;               /* LAMMPS handle */
    int         initialized;        /* mesh+solver initialized? */
} CouplingState;

/* ---------------------------------------------------------------------------
 * Utility macros
 * -------------------------------------------------------------------------*/
#define COORD(mesh, i, d)  ((mesh)->coord[3*(i) + (d)])
#define ETPL(mesh, e, n)   ((mesh)->etpl[8*(e) + (n)])

/* Safe memory allocation */
static inline void *mpm_calloc(size_t count, size_t size) {
    void *p = calloc(count, size);
    if (!p && count > 0) {
        fprintf(stderr, "mpm_calloc: out of memory (%zu × %zu)\n", count, size);
        exit(1);
    }
    return p;
}

static inline void *mpm_realloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size > 0) {
        fprintf(stderr, "mpm_realloc: out of memory (%zu)\n", size);
        exit(1);
    }
    return p;
}

#endif /* MPM_TYPES_H */
