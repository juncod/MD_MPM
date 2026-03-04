# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**AMPLE (A Material Point Learning Environment)** — MATLAB-based Material Point Method (MPM) framework for solving PDEs. Currently contains:
- `MPM_EXAMPLE_HEAT/` — mechanical heat conduction (baseline reference implementation)
- `MPM_EXAMPLE_ELECTRIC/` — 2D steady-state electric conduction (adapted from HEAT)
- `MPM_EXAMPLE_ELECTRIC_3D/` — 3D steady-state electric conduction (extended from 2D)
- `LAMMPS_EXAMPLE/` — Cu grain boundary MD simulations (independent, not integrated with MPM)

## Running Examples

```bash
# Run MATLAB examples (headless, no GUI)
matlab -batch "cd('E:/Project_MD_MPM/MPM_EXAMPLE_ELECTRIC'); ample_electric"
matlab -batch "cd('E:/Project_MD_MPM/MPM_EXAMPLE_ELECTRIC_3D'); ample_electric_3d"
matlab -batch "cd('E:/Project_MD_MPM/MPM_EXAMPLE_HEAT'); ample_heat"

# Visualize VTK output with Python
cd MPM_EXAMPLE_ELECTRIC/output && python visualize_vtk.py
cd MPM_EXAMPLE_ELECTRIC_3D/output && python visualize_vtk_3d.py
```

**Prerequisites**: `output/` directory must exist before running (VTK files are written there).

## Code Architecture

### Computational Flow (all examples follow this pattern)

```
ample_*.m (main driver)
  └── setupGrid_*.m         — define mesh, MPs, BCs, material constants
  └── [loadstep loop]
        └── elemMPinfo_*.m  — compute MP↔node connectivity + shape functions
        └── detExtForce_*.m — assemble external force/current vector
        └── [Newton-Raphson]
              └── detFDoFs_*.m    — identify free DOFs (non-Dirichlet)
              └── linSolve_*.m    — solve K·u = f with Dirichlet BCs
              └── detMPs_*.m      — assemble global K & internal forces
        └── updateMPs_*.m   — update material point state
        └── postPro_*.m     — write VTK output
```

### Core Data Structures

**`mpData` struct array** — one entry per material point:
- `.mpC` — coordinates (1×nD)
- `.vp`, `.vp0` — current/initial volume
- `.lp`, `.lp0` — GIMP domain half-lengths (nD×1)
- `.nIN` — connected node indices
- `.Svp`, `.dSvp` — basis values and gradients
- `.mCst` — material constants `[σ, β]` (conductivity, potential-dependence)
- `.phi`, `.gradPhi`, `.J` — electric potential, gradient, current density
- `.fp` — volumetric current source (external)

**`mesh` struct**:
- `.etpl` — element topology (nels × nen)
- `.coord` — node coordinates (nnodes × nD)
- `.bc` — Dirichlet BCs (nbcs × 2): `[node_index, prescribed_value]`
- `.h` — element sizes (nD × 1)

### Physics: Electric Conduction

Solves weak form of: **−∇·(σ∇φ) = 0**

Key assembly in `detMPs_Electric.m`:
```matlab
G = dNx;                    % ∇N  (gradient of basis)
D = sigma_c * eye(nD);      % conductivity tensor
Ke += (G' * D * G) * vol;   % conductance: ∫ ∇N'·σ·∇N dΩ
fint += (G' * D * G * phi_nodes) * vol;  % internal flux
J = -D * G * phi_nodes;     % Ohm's law: J = -σ∇φ
```

### Shape Functions

- **MPM** (`SvpMPM.m`): piecewise-linear basis per element
- **GIMP** (`SvpGIMP.m`): generalized interpolation — smoother, reduces noise at cell boundaries (controlled by `mpType=2` in setup)

`MPMbasis.m` wraps the 1D functions and applies them per dimension via tensor product.

### Adding a New Physics Example

Follow the pattern from ELECTRIC → ELECTRIC_3D:
1. Create folder with `ample_*.m`, `setup/`, `functions/`, `plotting/`, `output/`
2. Adapt `setupGrid_*.m` for domain, mesh, MPs, BCs, material constants
3. Create `*_Physics` variants of: `elemMPinfo`, `detExtForce`, `detFDoFs`, `linSolve`, `detMPs`, `updateMPs`
4. Create `makeVtkMP_*.m` writer and `postPro_*.m` orchestrator
5. Can reuse shared utilities: `formCoord2D/3D.m`, `MPMbasis.m`, `SvpMPM/GIMP.m`, `nodesForMP.m`, `elemForMP.m`

## VTK Output Format

ASCII VTK Unstructured Grid — written by `makeVtkMP_Electric.m`:
- `mesh.vtk` — background FE mesh (QUAD or HEX cells)
- `mpData_i.vtk` — material points at loadstep `i` with `SCALARS phi` and `VECTORS J`

## Problem Parameters (Electric 2D)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 42×40 | `nelsx×nelsy` elements (lx=42, ly=40) |
| Conductivity σ | 1.0 | Uniform |
| MP type | GIMP (2) | 2×2 MPs per element |
| BC left | φ=0 | Ground at x=0 |
| BC right | φ=1 | 1V at x=42 |
| Hole | center=[21,20], r=5 | MPs inside removed |
| Expected J_mean | ~1/42 ≈ 0.02381 | Theoretical value |
| Current crowding | ~2× at hole flanks | Theoretical max |
