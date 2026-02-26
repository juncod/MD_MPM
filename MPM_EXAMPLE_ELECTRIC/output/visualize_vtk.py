"""
Visualize MPM Electric Field VTK output using matplotlib.
Reads mpData_1.vtk (material points) and mesh.vtk (background mesh).
"""
import matplotlib
matplotlib.use('Agg')   # non-interactive: save to file without GUI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os

# ── VTK reader ─────────────────────────────────────────────────────────────
def read_mp_vtk(path):
    """Parse ASCII VTK unstructured grid written by makeVtkMP_Electric."""
    with open(path, 'r') as f:
        lines = f.readlines()

    pts, phi, J = [], [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('POINTS'):
            n = int(line.split()[1])
            for _ in range(n):
                i += 1
                pts.append(list(map(float, lines[i].split())))
        elif line.startswith('SCALARS phi'):
            i += 1          # skip LOOKUP_TABLE line
            while i + 1 < len(lines) and not lines[i+1].strip().startswith('VECTORS'):
                i += 1
                v = lines[i].strip()
                if v:
                    phi.append(float(v))
        elif line.startswith('VECTORS J'):
            while i + 1 < len(lines):
                i += 1
                v = lines[i].strip()
                if not v:
                    break
                J.append(list(map(float, v.split())))
        i += 1

    pts = np.array(pts)
    phi = np.array(phi)
    J   = np.array(J)
    return pts[:, 0], pts[:, 1], phi, J[:, 0], J[:, 1]


def read_mesh_vtk(path):
    """Parse ASCII VTK mesh (background grid nodes + cells)."""
    with open(path, 'r') as f:
        lines = f.readlines()

    nodes, cells = [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('POINTS'):
            n = int(line.split()[1])
            for _ in range(n):
                i += 1
                nodes.append(list(map(float, lines[i].split()))[:2])
        elif line.startswith('CELLS'):
            n = int(line.split()[1])
            for _ in range(n):
                i += 1
                vals = list(map(int, lines[i].split()))
                cells.append(vals[1:vals[0]+1])   # skip the count prefix
        i += 1
    return np.array(nodes), cells


# ── Load data ──────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
mp_vtk   = os.path.join(script_dir, 'mpData_1.vtk')
mesh_vtk = os.path.join(script_dir, 'mesh.vtk')

x, y, phi, Jx, Jy = read_mp_vtk(mp_vtk)
mesh_nodes, mesh_cells = read_mesh_vtk(mesh_vtk)

J_mag = np.sqrt(Jx**2 + Jy**2)
J0    = 1 / 42          # theoretical uniform current density

print(f"Material points: {len(x)}")
print(f"phi  : min={phi.min():.4f}  max={phi.max():.4f}")
print(f"|J|  : min={J_mag.min():.5f}  max={J_mag.max():.5f}  mean={J_mag.mean():.5f}")
print(f"Max/J0 (crowding ratio) = {J_mag.max()/J0:.2f}")

# ── Helper: draw mesh & hole ───────────────────────────────────────────────
def draw_mesh(ax, alpha=0.15):
    for cell in mesh_cells:
        poly = mesh_nodes[cell]
        poly = np.vstack([poly, poly[0]])
        ax.plot(poly[:, 0], poly[:, 1], 'k-', lw=0.3, alpha=alpha)

def draw_hole(ax, color='white', lw=1.5):
    hole = Circle((21, 20), 5, fill=True, facecolor=color,
                  edgecolor='red', linewidth=lw, zorder=5)
    ax.add_patch(hole)

# ── Figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('MPM Electric Field — 42×40 Conductor Plate with Circular Hole (r=5)',
             fontsize=13, fontweight='bold')

mp_size = 4   # scatter marker size

# ── Panel 1: Electric potential φ ─────────────────────────────────────────
ax = axes[0]
sc = ax.scatter(x, y, c=phi, s=mp_size, cmap='RdYlBu_r',
                vmin=0, vmax=1, rasterized=True)
draw_mesh(ax)
draw_hole(ax)
cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
cb.set_label('φ (V)', fontsize=11)
ax.set_title('Electric Potential φ', fontsize=12)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_aspect('equal'); ax.set_xlim(-0.5, 42.5); ax.set_ylim(-0.5, 40.5)
# Mark electrodes
ax.axvline(0,  color='blue',  lw=2, ls='--', label='φ=0 (ground)', alpha=0.7)
ax.axvline(42, color='red',   lw=2, ls='--', label='φ=1V',         alpha=0.7)
ax.legend(fontsize=8, loc='upper left')

# ── Panel 2: Current density magnitude |J| ────────────────────────────────
ax = axes[1]
sc2 = ax.scatter(x, y, c=J_mag, s=mp_size, cmap='hot_r',
                 vmin=0, vmax=J_mag.max(), rasterized=True)
draw_mesh(ax)
draw_hole(ax)
cb2 = plt.colorbar(sc2, ax=ax, fraction=0.046, pad=0.04)
cb2.set_label('|J| (A/m²)', fontsize=11)
ax.set_title('Current Density Magnitude |J|', fontsize=12)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_aspect('equal'); ax.set_xlim(-0.5, 42.5); ax.set_ylim(-0.5, 40.5)
# Annotate crowding spots
ax.annotate(f'max={J_mag.max():.4f}\n({J_mag.max()/J0:.1f}×J₀)',
            xy=(21, 20+5.5), xytext=(28, 28),
            arrowprops=dict(arrowstyle='->', color='cyan', lw=1.5),
            color='cyan', fontsize=9, fontweight='bold')

# ── Panel 3: Current vector field (quiver) near hole ─────────────────────
ax = axes[2]
# Full-field scatter (faded)
ax.scatter(x, y, c=J_mag, s=mp_size, cmap='hot_r',
           vmin=0, vmax=J_mag.max(), alpha=0.4, rasterized=True)
draw_mesh(ax)
draw_hole(ax, color='lightgray')

# Quiver: subsample for clarity
step = max(1, len(x) // 500)
qx, qy = x[::step], y[::step]
qJx, qJy = Jx[::step], Jy[::step]
q = ax.quiver(qx, qy, qJx, qJy, J_mag[::step],
              cmap='plasma', scale=0.8, scale_units='xy',
              width=0.004, headwidth=3, headlength=4, zorder=4)
plt.colorbar(q, ax=ax, fraction=0.046, pad=0.04, label='|J| (A/m²)')
ax.set_title('Current Density Vectors J\n(current crowding visible near hole)', fontsize=11)
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_aspect('equal'); ax.set_xlim(-0.5, 42.5); ax.set_ylim(-0.5, 40.5)

# ── Save ───────────────────────────────────────────────────────────────────
out_png = os.path.join(script_dir, 'electric_result.png')
plt.tight_layout()
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out_png}")

# ── Extra: zoom view around hole ───────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 6))
fig2.suptitle('Zoom: Current Crowding Around Circular Hole', fontsize=13, fontweight='bold')

# Zoom region
xz1, xz2, yz1, yz2 = 10, 32, 10, 32
mask = (x > xz1) & (x < xz2) & (y > yz1) & (y < yz2)

ax = axes2[0]
sc = ax.scatter(x[mask], y[mask], c=phi[mask], s=12, cmap='RdYlBu_r',
                vmin=0, vmax=1)
draw_hole(ax, color='white')
cb = plt.colorbar(sc, ax=ax); cb.set_label('φ (V)')
ax.set_title('Potential φ (zoom)', fontsize=12)
ax.set_xlim(xz1, xz2); ax.set_ylim(yz1, yz2); ax.set_aspect('equal')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

ax = axes2[1]
sc2 = ax.scatter(x[mask], y[mask], c=J_mag[mask], s=12, cmap='hot_r',
                 vmin=0, vmax=J_mag.max())
draw_hole(ax, color='lightgray')
# Dense quiver in zoom area
step_z = max(1, mask.sum() // 200)
idx = np.where(mask)[0][::step_z]
ax.quiver(x[idx], y[idx], Jx[idx], Jy[idx], J_mag[idx],
          cmap='winter', scale=0.6, scale_units='xy',
          width=0.006, headwidth=3, headlength=4, zorder=5)
cb2 = plt.colorbar(sc2, ax=ax); cb2.set_label('|J| (A/m²)')
ax.set_title('Current Density |J| + vectors (zoom)', fontsize=12)
ax.set_xlim(xz1, xz2); ax.set_ylim(yz1, yz2); ax.set_aspect('equal')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

out_zoom = os.path.join(script_dir, 'electric_result_zoom.png')
plt.tight_layout()
plt.savefig(out_zoom, dpi=150, bbox_inches='tight')
print(f"Saved: {out_zoom}")

# plt.show()  # disabled for non-interactive batch run
