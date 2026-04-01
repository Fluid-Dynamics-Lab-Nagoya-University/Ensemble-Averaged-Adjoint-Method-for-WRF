"""
Toy model for illustrating the proposed method.

by Shan Jiang, FDL, Nagoya University
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import minimum_filter
import matplotlib as mpl
import os

# =============================================================================
# ── PARAMETRIC CONTROL PANEL ─────────────────────────────────────────────────
# =============================================================================
n_ensemble_list = [100, 200, 300, 400]  # list of ensemble sizes to run

global_scale = 1.2          # global scale factor: font size, marker size, arrow size/linewidth/length

tri_display_radius = 1      # radius around the initial point within which local minima markers are shown
base_font_size = 16 * global_scale
selected_cmap = plt.cm.coolwarm_r  # colormap for the landscape surface (reversed cool-warm)

ideal_arrow_lw = 2.0 * global_scale;  ideal_arrow_ms = 20 * global_scale
avg_arrow_lw = 3.5 * global_scale;    avg_arrow_ms = 20 * global_scale
det_arrow_lw = 3.0 * global_scale;    det_arrow_ms = 20 * global_scale

member_arrow_color = 'white'
member_arrow_alpha = 0.7
member_arrow_lw = 1.0 * global_scale
member_arrow_ms = 14 * global_scale
member_point_ms = 4 * global_scale
member_arrow_scale = 0.1 * global_scale

x0_pos = (-0.3, -0.3)
global_min_target = (1.8, 1.8)
# =============================================================================

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = base_font_size

if not os.path.exists('./figs'):
    os.makedirs('./figs')

DPI = 150
CBAR_LABEL = r'$\mathcal{J}$'

NPTS = 500
x = np.linspace(-3, 3, NPTS)
y = np.linspace(-3, 3, NPTS)
X, Y = np.meshgrid(x, y)

def make_grf(shape, seed, power_exp, amplitude):
    rng = np.random.default_rng(seed)
    F = np.fft.fft2(rng.standard_normal(shape))
    ky = np.fft.fftfreq(shape[0])
    kx = np.fft.fftfreq(shape[1])
    KX, KY = np.meshgrid(kx, ky)
    K = np.sqrt(KX**2 + KY**2)
    K[0, 0] = 1.0
    power = K ** power_exp
    power[0, 0] = 0.0
    field = np.real(np.fft.ifft2(F * np.sqrt(power)))
    return amplitude * field / field.std()

Z = 0.02 * (X**2 + Y**2)
Z -= 40.0 * np.exp(-((X - global_min_target[0])**2
                      + (Y - global_min_target[1])**2) / 4.0)
Z -= 3.0 * np.exp(-((X - (-0.42))**2 + (Y - (-0.42))**2) / 0.03)
for cx_, cy_, d, w in [
    (-2.0,  1.0, 1.5, 0.30), ( 1.0, -2.0, 1.5, 0.30),
    (-2.0, -1.5, 1.2, 0.25), (-1.5,  2.0, 1.2, 0.25),
    ( 2.3, -0.8, 1.0, 0.25), ( 0.0,  2.5, 0.8, 0.20),
]:
    Z -= d * np.exp(-((X - cx_)**2 + (Y - cy_)**2) / w)

GRF_SEED = 633
Z += make_grf(X.shape, GRF_SEED,       -3.0, 6.0)
Z += make_grf(X.shape, GRF_SEED+1000,  -2.5, 3.0)
Z += make_grf(X.shape, GRF_SEED+2000,  -2.0, 1.2)

min_idx = np.argmin(Z)
iy_min, ix_min = np.unravel_index(min_idx, Z.shape)
true_global_min_x = X[iy_min, ix_min]
true_global_min_y = Y[iy_min, ix_min]

neighborhood = 20
local_min_mask = (Z == minimum_filter(Z, size=neighborhood))
iy_all, ix_all = np.where(local_min_mask)
true_local_min_pts = []
for iy, ix in zip(iy_all, ix_all):
    px, py = X[iy, ix], Y[iy, ix]
    dist_to_global = np.sqrt((px - true_global_min_x)**2 + (py - true_global_min_y)**2)
    dist_to_start = np.sqrt((px - x0_pos[0])**2 + (py - x0_pos[1])**2)
    if dist_to_global > 0.6 and dist_to_start <= tri_display_radius:
        true_local_min_pts.append((px, py))

def grad_at(Z_mat, px, py):
    ix_p = np.argmin(np.abs(x - px))
    iy_p = np.argmin(np.abs(y - py))
    ix_p = np.clip(ix_p, 1, len(x) - 2)
    iy_p = np.clip(iy_p, 1, len(y) - 2)
    dx = (Z_mat[iy_p, ix_p + 1] - Z_mat[iy_p, ix_p - 1]) / (2 * (x[1] - x[0]))
    dy = (Z_mat[iy_p + 1, ix_p] - Z_mat[iy_p - 1, ix_p]) / (2 * (y[1] - y[0]))
    return -dx, -dy

vmin, vmax = Z.min(), Z.max() * 0.5
arrow_len = 1.05 * global_scale

def rot90cw(px, py):
    return py, -px

ideal_dir = np.array([global_min_target[0] - x0_pos[0],
                      global_min_target[1] - x0_pos[1]])
ideal_dir = ideal_dir / np.linalg.norm(ideal_dir)

sigma_pert_list = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

# =============================================================================
# ── FIGURE GENERATION + DIAGNOSTIC (merged) ──────────────────────────────────
# =============================================================================
X_r = Y
Y_r = -X

# ── precompute deterministic arrow direction (after rotation) ──
gx_det, gy_det = grad_at(Z, x0_pos[0], x0_pos[1])
nm_det = np.sqrt(gx_det**2 + gy_det**2)
det_gxr, det_gyr = rot90cw(gx_det / nm_det * arrow_len, gy_det / nm_det * arrow_len)

det_dir = np.array([gx_det, gy_det])
det_dir_n = det_dir / np.linalg.norm(det_dir)
det_angle = np.degrees(np.arccos(np.clip(np.dot(det_dir_n, ideal_dir), -1, 1)))

diag_path = './figs/toy_diagnostic.txt'
with open(diag_path, 'w') as f_diag:
    def log(msg=''):
        print(msg)
        f_diag.write(msg + '\n')

    log("=" * 60)
    log("Diagnostic: angle between ensemble-averaged gradient")
    log("            and ideal direction toward global minimum")
    log("=" * 60)
    log()
    log(f"  x0           = {x0_pos}")
    log(f"  global_min   = ({true_global_min_x:.2f}, {true_global_min_y:.2f})")
    log(f"  local mins displayed = {len(true_local_min_pts)}")
    log()
    log(f"  Deterministic  angle_to_ideal = {det_angle:6.1f} deg")
    log()

    for n_ensemble in n_ensemble_list:
        log(f"--- N_ensemble = {n_ensemble} ---")

        for sp in sigma_pert_list:
            label = (f'($N$={n_ensemble}, $\\sigma_{{\\mathrm{{pert}}}}$={sp})')
            fname = f'landscape_b_N{n_ensemble}_S{sp}.png'

            rng_panel = np.random.default_rng(42)

            fig, ax = plt.subplots(figsize=(6.2, 5.6))
            cf = ax.contourf(X_r, Y_r, Z, levels=65, cmap=selected_cmap,
                             vmin=vmin, vmax=vmax, alpha=0.95)
            ax.contour(X_r, Y_r, Z, levels=40, colors='white',
                       linewidths=0.25*global_scale, alpha=0.25)

            for tx, ty in true_local_min_pts:
                tx_r, ty_r = rot90cw(tx, ty)
                ax.plot(tx_r, ty_r, 'w^', ms=5*global_scale, zorder=8,
                        markeredgecolor='0.3', alpha=0.8)

            gmin_xr, gmin_yr = rot90cw(true_global_min_x, true_global_min_y)
            ax.plot(gmin_xr, gmin_yr, '*', color='gold',
                    ms=16*global_scale, zorder=20, markeredgecolor='black', markeredgewidth=0.8*global_scale)

            x0r, y0r = rot90cw(x0_pos[0], x0_pos[1])
            ax.plot(x0r, y0r, 'ko', ms=9*global_scale, zorder=22)
            ax.text(x0r - 0.12, y0r - 0.3, r'$\mathbf{x}^0$',
                    fontweight='bold', zorder=25)

            tx_ideal = true_global_min_x - x0_pos[0]
            ty_ideal = true_global_min_y - x0_pos[1]
            tn = np.sqrt(tx_ideal**2 + ty_ideal**2)
            id_dx_r, id_dy_r = rot90cw(tx_ideal / tn * arrow_len,
                                       ty_ideal / tn * arrow_len)

            # ── ideal direction (dashed green) ──
            ax.annotate('',
                        xy=(x0r + id_dx_r, y0r + id_dy_r),
                        xytext=(x0r, y0r),
                        arrowprops=dict(arrowstyle='->', color='#00CC77',
                                        lw=ideal_arrow_lw, ls='--', alpha=0.5,
                                        mutation_scale=ideal_arrow_ms),
                        zorder=28)

            # ── deterministic adjoint (red arrow) ──
            ax.annotate('',
                        xy=(x0r + det_gxr, y0r + det_gyr),
                        xytext=(x0r, y0r),
                        arrowprops=dict(arrowstyle='->', color='#FF3333',
                                        lw=det_arrow_lw,
                                        mutation_scale=det_arrow_ms),
                        zorder=29)

            # ── ensemble members ──
            ensemble_grads = []
            for i in range(n_ensemble):
                ex = x0_pos[0] + rng_panel.standard_normal() * sp
                ey = x0_pos[1] + rng_panel.standard_normal() * sp
                egx, egy = grad_at(Z, ex, ey)
                ensemble_grads.append((egx, egy))

                mag = np.sqrt(egx**2 + egy**2)
                if mag > 0:
                    log_len = np.log1p(mag) * member_arrow_scale
                    dir_xr, dir_yr = rot90cw(egx / mag * log_len,
                                             egy / mag * log_len)
                else:
                    dir_xr, dir_yr = 0, 0

                exr, eyr = rot90cw(ex, ey)
                ax.plot(exr, eyr, 'o', color=member_arrow_color,
                        ms=member_point_ms, alpha=member_arrow_alpha, zorder=10)
                ax.annotate('',
                            xy=(exr + dir_xr, eyr + dir_yr),
                            xytext=(exr, eyr),
                            arrowprops=dict(arrowstyle='->',
                                            color=member_arrow_color,
                                            alpha=member_arrow_alpha,
                                            lw=member_arrow_lw,
                                            mutation_scale=member_arrow_ms),
                            zorder=11)

            # ── ensemble-averaged direction (solid green) ──
            avg_gx = np.mean([g[0] for g in ensemble_grads])
            avg_gy = np.mean([g[1] for g in ensemble_grads])
            avg_nm = np.sqrt(avg_gx**2 + avg_gy**2)
            final_gxr, final_gyr = rot90cw(avg_gx / avg_nm * arrow_len,
                                           avg_gy / avg_nm * arrow_len)

            ax.annotate('',
                        xy=(x0r + final_gxr, y0r + final_gyr),
                        xytext=(x0r, y0r),
                        arrowprops=dict(arrowstyle='->', color='#00CC77',
                                        lw=avg_arrow_lw,
                                        mutation_scale=avg_arrow_ms),
                        zorder=30)

            # ── diagnostic for this (N, sigma) pair ──
            avg_dir = np.array([avg_gx, avg_gy])
            avg_dir_norm = avg_dir / np.linalg.norm(avg_dir)
            angle_to_ideal = np.degrees(np.arccos(np.clip(np.dot(avg_dir_norm, ideal_dir), -1, 1)))
            angle_to_det = np.degrees(np.arccos(np.clip(np.dot(avg_dir_norm, det_dir_n), -1, 1)))
            log(f"  sigma_pert={sp:<6.3f}  angle_to_ideal={angle_to_ideal:6.1f} deg  angle_to_det={angle_to_det:5.1f} deg")

            ax.set_xlim(-3, 3)
            ax.set_ylim(-3, 3)
            ax.set_aspect('equal')
            ax.set_title(label, fontweight='bold', pad=12)
            ax.tick_params(left=False, bottom=False,
                           labelleft=False, labelbottom=False)

            cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(CBAR_LABEL)
            cbar.set_ticks([])

            plt.tight_layout()
            out_fname = f'./figs/N{n_ensemble}_{fname.replace(".png", f"_dpi{DPI}.png")}'
            plt.savefig(out_fname, dpi=DPI, bbox_inches='tight')
            plt.close()
            print(f'Saved: {out_fname}')

        log()

    log("=" * 60)

print(f"Saved: {diag_path}")