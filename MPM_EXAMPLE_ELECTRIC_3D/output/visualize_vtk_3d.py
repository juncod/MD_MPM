"""
Visualize 3D MPM Electric Field VTK output using matplotlib.
Reads mpData_1.vtk (material points) and produces 2D slice plots.
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os

# ── VTK reader ─────────────────────────────────────────────────────────────
def read_mp_vtk_3d(path):
    """Parse ASCII VTK unstructured grid (3D) from makeVtkMP_Electric_3d."""
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
            i += 1  # skip LOOKUP_TABLE line
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
    return pts, phi, J


# ── Load data ───────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
pts, phi, J = read_mp_vtk_3d(os.path.join(script_dir, 'mpData_1.vtk'))
x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
Jx, Jy, Jz = J[:, 0], J[:, 1], J[:, 2]
J_mag = np.sqrt(Jx**2 + Jy**2 + Jz**2)

J0 = 1.0 / 20.0   # theoretical uniform current density (1V / 20m)

print(f"Material points: {len(x)}")
print(f"phi  : min={phi.min():.5f}  max={phi.max():.5f}")
print(f"|J|  : min={J_mag.min():.6f}  max={J_mag.max():.6f}  mean={J_mag.mean():.6f}")
print(f"J0   : {J0:.6f}")
print(f"Max/J0 (crowding ratio) = {J_mag.max()/J0:.3f}  (theory 3D sphere: 1.5)")

# ── Slice selection (mid-plane through sphere center) ───────────────────────
mid = np.array([10.0, 10.0, 10.0])
tol = 0.6   # half-thickness of slice

# XY-plane slice (z ≈ 10)
mXY = np.abs(z - mid[2]) < tol
# XZ-plane slice (y ≈ 10)
mXZ = np.abs(y - mid[1]) < tol
# YZ-plane slice (x ≈ 10)
mYZ = np.abs(x - mid[0]) < tol

def draw_sphere_circle(ax, plane='xy', r=3.0, center=(10,10)):
    """Draw sphere cross-section circle on a 2D slice."""
    c = Circle(center, r, fill=False, edgecolor='cyan', linewidth=1.5,
               linestyle='--', zorder=5)
    ax.add_patch(c)


# ── Figure 1: 3-panel slice overview ────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('MPM 3D Electric Field — 20×20×20 Domain with Spherical Hole (r=3)\n'
             'Slice plots through sphere center', fontsize=13, fontweight='bold')

ms = 6  # marker size

# --- Row 0: Electric potential φ ---
for col, (mask, sl_x, sl_y, lx, ly, xlab, ylab, title, cx, cy) in enumerate([
    (mXY, x, y, (0,20), (0,20), 'x', 'y', 'φ — XY slice (z=10)', 10, 10),
    (mXZ, x, z, (0,20), (0,20), 'x', 'z', 'φ — XZ slice (y=10)', 10, 10),
    (mYZ, y, z, (0,20), (0,20), 'y', 'z', 'φ — YZ slice (x=10)', 10, 10),
]):
    ax = axes[0, col]
    sc = ax.scatter(sl_x[mask], sl_y[mask], c=phi[mask], s=ms,
                    cmap='RdYlBu_r', vmin=0, vmax=1, rasterized=True)
    draw_sphere_circle(ax, r=3.0, center=(cx, cy))
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label('φ (V)')
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlab + ' (m)'); ax.set_ylabel(ylab + ' (m)')
    ax.set_xlim(*lx); ax.set_ylim(*ly); ax.set_aspect('equal')
    if col == 0:
        ax.axvline(0,  color='blue', lw=1.5, ls='--', alpha=0.7, label='φ=0')
        ax.axvline(20, color='red',  lw=1.5, ls='--', alpha=0.7, label='φ=1V')
        ax.legend(fontsize=8)

# --- Row 1: Current density |J| ---
vmax_J = J_mag.max()
for col, (mask, sl_x, sl_y, lx, ly, xlab, ylab, title, cx, cy, qx2, qy2, qJx2, qJy2) in enumerate([
    (mXY, x, y, (0,20), (0,20), 'x', 'y', '|J| — XY slice (z=10)', 10, 10,
     x[mXY], y[mXY], Jx[mXY], Jy[mXY]),
    (mXZ, x, z, (0,20), (0,20), 'x', 'z', '|J| — XZ slice (y=10)', 10, 10,
     x[mXZ], z[mXZ], Jx[mXZ], Jz[mXZ]),
    (mYZ, y, z, (0,20), (0,20), 'y', 'z', '|J| — YZ slice (x=10)', 10, 10,
     y[mYZ], z[mYZ], Jy[mYZ], Jz[mYZ]),
]):
    ax = axes[1, col]
    sc = ax.scatter(sl_x[mask], sl_y[mask], c=J_mag[mask], s=ms,
                    cmap='hot_r', vmin=0, vmax=vmax_J, rasterized=True)
    draw_sphere_circle(ax, r=3.0, center=(cx, cy))
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label('|J| (A/m²)')
    # Quiver subsample
    step = max(1, len(qx2) // 300)
    ax.quiver(qx2[::step], qy2[::step], qJx2[::step], qJy2[::step],
              scale=1.0, scale_units='xy', width=0.003,
              headwidth=3, headlength=4, color='gray', alpha=0.6, zorder=4)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlab + ' (m)'); ax.set_ylabel(ylab + ' (m)')
    ax.set_xlim(*lx); ax.set_ylim(*ly); ax.set_aspect('equal')
    ax.annotate(f'max={J_mag.max():.4f}\n({J_mag.max()/J0:.2f}×J₀)',
                xy=(cx, cy + 3.5), xytext=(cx + 4, cy + 7),
                arrowprops=dict(arrowstyle='->', color='cyan', lw=1.5),
                color='cyan', fontsize=9, fontweight='bold')

out_png = os.path.join(script_dir, 'electric_result_3d.png')
plt.tight_layout()
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out_png}")

# ── Figure 2: Zoom around sphere (XY slice) ─────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6))
fig2.suptitle('Zoom: Current Crowding Around Spherical Hole (XY slice, z≈10)',
              fontsize=13, fontweight='bold')

xz1, xz2, yz1, yz2 = 5, 15, 5, 15
mzoom = mXY & (x > xz1) & (x < xz2) & (y > yz1) & (y < yz2)

# Panel 1: phi zoom
ax = axes2[0]
sc = ax.scatter(x[mzoom], y[mzoom], c=phi[mzoom], s=20, cmap='RdYlBu_r', vmin=0, vmax=1)
draw_sphere_circle(ax, r=3.0, center=(10, 10))
plt.colorbar(sc, ax=ax).set_label('φ (V)')
ax.set_title('Potential φ (zoom, XY slice)', fontsize=12)
ax.set_xlim(xz1, xz2); ax.set_ylim(yz1, yz2); ax.set_aspect('equal')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

# Panel 2: |J| zoom
ax = axes2[1]
sc2 = ax.scatter(x[mzoom], y[mzoom], c=J_mag[mzoom], s=20, cmap='hot_r', vmin=0, vmax=vmax_J)
draw_sphere_circle(ax, r=3.0, center=(10, 10))
plt.colorbar(sc2, ax=ax).set_label('|J| (A/m²)')
ax.set_title('Current Density |J| (zoom)', fontsize=12)
ax.set_xlim(xz1, xz2); ax.set_ylim(yz1, yz2); ax.set_aspect('equal')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

# Panel 3: J vectors zoom (dense quiver)
ax = axes2[2]
sc3 = ax.scatter(x[mzoom], y[mzoom], c=J_mag[mzoom], s=20, cmap='hot_r',
                 vmin=0, vmax=vmax_J, alpha=0.5)
draw_sphere_circle(ax, r=3.0, center=(10, 10))
step_z = max(1, mzoom.sum() // 200)
idx = np.where(mzoom)[0][::step_z]
ax.quiver(x[idx], y[idx], Jx[idx], Jy[idx], J_mag[idx],
          cmap='winter', scale=0.5, scale_units='xy',
          width=0.008, headwidth=3, headlength=4, zorder=5)
plt.colorbar(sc3, ax=ax).set_label('|J| (A/m²)')
ax.set_title('J vectors — crowding near sphere equator', fontsize=12)
ax.set_xlim(xz1, xz2); ax.set_ylim(yz1, yz2); ax.set_aspect('equal')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')

out_zoom = os.path.join(script_dir, 'electric_result_3d_zoom.png')
plt.tight_layout()
plt.savefig(out_zoom, dpi=150, bbox_inches='tight')
print(f"Saved: {out_zoom}")

# ── Figure 3: J_mag profile along x-axis (y=z=10) ──────────────────────────
fig3, ax3 = plt.subplots(figsize=(10, 5))
# MPs along the flow axis (x-direction, y≈10, z≈10)
m_axis = (np.abs(y - 10) < 0.3) & (np.abs(z - 10) < 0.3)
xs_axis = x[m_axis]
Jm_axis = J_mag[m_axis]
idx_sort = np.argsort(xs_axis)
ax3.plot(xs_axis[idx_sort], Jm_axis[idx_sort], 'b-o', ms=3, lw=1.5, label='|J| on axis (y=z=10)')
ax3.axhline(J0, color='gray', lw=1.5, ls='--', label=f'J₀ = {J0:.4f} (uniform)')
# Mark sphere boundary on axis
ax3.axvspan(7, 13, alpha=0.15, color='red', label='Sphere region (void)')
ax3.set_xlabel('x (m)', fontsize=12)
ax3.set_ylabel('|J| (A/m²)', fontsize=12)
ax3.set_title('Current Density Along x-axis (y=z=10)\n— sphere hole creates current shadow', fontsize=12)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.4)

out_prof = os.path.join(script_dir, 'electric_result_3d_profile.png')
plt.tight_layout()
plt.savefig(out_prof, dpi=150, bbox_inches='tight')
print(f"Saved: {out_prof}")
