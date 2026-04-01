"""
Plots 3-D toy model for illustrating the proposed method.

by Shan Jiang, FDL, Nagoya University
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import os

# =============================================================================
# ── PARAMETRIC CONTROL PANEL ─────────────────────────────────────────────────
# =============================================================================
ELEV        = 30           # 3D view elevation angle (degrees)
AZIM        = 390          # 3D view azimuth angle (degrees)
CMAP        = 'coolwarm_r' # colormap for the landscape surface
FONT_SIZE   = 16           # global font size
FONT_FAMILY = 'Helvetica'

FIG_W       = 6.5          # figure width (inches)
FIG_H       = 5.5          # figure height (inches)

# position and size of the 3D axes within the figure
SUBPLOT_RECT = [0.001, 0.001, 0.999, 0.999]
# =============================================================================

RSTRIDE     = 2
CSTRIDE     = 2
ALPHA_3D    = 0.9

X0_COLOR    = 'black'
XMIN_COLOR  = 'gold'
X0_SIZE     = 8
XMIN_SIZE   = 14

ARROW_COLOR = 'black'
ARROW_LW    = 2
ARROW_MS    = 20
ARROW_SHRINK = 0.08
ARROW_DASH  = (3.2, 3.2)

# =============================================================================

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = [FONT_FAMILY, 'Helvetica']
mpl.rcParams['font.size'] = FONT_SIZE

os.makedirs('./figs', exist_ok=True)

x0_pos = (-0.3, -0.3)
global_min_target = (1.8, 1.8)
GRF_SEED = 633
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
Z -= 40.0 * np.exp(-((X - global_min_target[0])**2 + (Y - global_min_target[1])**2) / 4.0)
Z -= 3.0 * np.exp(-((X - (-0.42))**2 + (Y - (-0.42))**2) / 0.03)
for cx_, cy_, d, w in [
    (-2.0,  1.0, 1.5, 0.30), ( 1.0, -2.0, 1.5, 0.30),
    (-2.0, -1.5, 1.2, 0.25), (-1.5,  2.0, 1.2, 0.25),
    ( 2.3, -0.8, 1.0, 0.25), ( 0.0,  2.5, 0.8, 0.20),
]:
    Z -= d * np.exp(-((X - cx_)**2 + (Y - cy_)**2) / w)
Z += make_grf(X.shape, GRF_SEED,      -3.0, 6.0)
Z += make_grf(X.shape, GRF_SEED+1000, -2.5, 3.0)
Z += make_grf(X.shape, GRF_SEED+2000, -2.0, 1.2)

min_idx = np.argmin(Z)
iy_min, ix_min = np.unravel_index(min_idx, Z.shape)
xmin_x = X[iy_min, ix_min]
xmin_y = Y[iy_min, ix_min]
xmin_z = Z[iy_min, ix_min]
x0_iz  = np.argmin(np.abs(y - x0_pos[1]))
x0_ix  = np.argmin(np.abs(x - x0_pos[0]))
x0_z   = Z[x0_iz, x0_ix]

# =============================================================================
# ── 3D LANDSCAPE FIGURE ───────────────────────────────────────────────────────
# =============================================================================
fig = plt.figure(figsize=(FIG_W, FIG_H))
ax1 = fig.add_axes(SUBPLOT_RECT, projection='3d')

ax1.plot_surface(X, Y, Z, cmap=CMAP, alpha=ALPHA_3D,
                 rstride=RSTRIDE, cstride=CSTRIDE,
                 linewidth=0, antialiased=False)

ax1.view_init(elev=ELEV, azim=AZIM)
ax1.set_zlim(Z.min(), Z.max())
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel(r'$\mathcal{J}$')

fig.canvas.draw()

def to_fig_coords(ax, x, y, z):
    x2d, y2d, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    screen = ax.transData.transform((x2d, y2d))
    return fig.transFigure.inverted().transform(screen)

p0f = to_fig_coords(ax1, x0_pos[0], x0_pos[1], x0_z)
p1f = to_fig_coords(ax1, xmin_x, xmin_y, xmin_z)

vec = p1f - p0f
p0f_s = p0f + vec * ARROW_SHRINK
p1f_s = p1f - vec * ARROW_SHRINK

ax_ov = fig.add_axes([0, 0, 1, 1], facecolor='none')
ax_ov.set_xlim(0, 1)
ax_ov.set_ylim(0, 1)
ax_ov.axis('off')

ax_ov.annotate('', xy=p1f_s, xytext=p0f_s,
               arrowprops=dict(arrowstyle='->', color=ARROW_COLOR,
                               lw=ARROW_LW, mutation_scale=ARROW_MS,
                               linestyle=(0, ARROW_DASH)))

ax_ov.plot(*p0f, 'o', color=X0_COLOR, markersize=X0_SIZE,
           markeredgecolor='black', markeredgewidth=0.8)
ax_ov.plot(*p1f, '*', color=XMIN_COLOR, markersize=XMIN_SIZE,
           markeredgecolor='black', markeredgewidth=0.8)

fig.savefig('./figs/landscape_3d.png', dpi=300, bbox_inches='tight')

plt.show()