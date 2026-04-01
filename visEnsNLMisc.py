#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mode 1: Read ensemble-averaged A_QVAPOR from .dat, normalize, compute
         ratio = (A_QVAPOR_norm * 0.02 * 0.5) / QVAPOR_sfc, and plot.
Mode 2: Read NL terminal RAINNC and AD t=6 RAINNC (100 members each),
         compute ensemble averages with 6 cores, plot NL_mean - AD_mean.
by Shan Jiang, FDL, Nagoya University
"""
import multiprocessing
multiprocessing.set_start_method('fork', force=True)
import numpy as np
from netCDF4 import Dataset
from multiprocessing import Pool
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from plotting_tools import draw_small_box

# ==========================
# Switch:  1 = ratio plot,  2 = RAINNC diff plot
# ==========================
mode = 2

# ==========================
# Parameters
# ==========================
ng   = 0.3   # background noise rate; 0 = single-member mode (reads only member 1)
ig   = 0.5
sg   = 0.3

TARGET_NORM = 0.095222  # target Frobenius norm for A_QVAPOR normalization


# suffix_switch: 1 = no suffix, 2 = append "N" after ig, 3 = append "P" after ig
suffix_switch = 1
_ig_suffix = {1: '', 2: 'N', 3: 'P'}.get(suffix_switch, '')




# ==========================
# Common paths
# ==========================
geo_file_path = "./clean/wrfinput_d01_clean"

# --- Mode 1 paths ---
dat_path = f"./dats/A_QVAPOR_mean_absG_i0M_NOV0.02_sg{sg}.dat"

# --- Mode 2 paths ---
nl_dir  = f"./wrfout/absG_i0M_NOV0.02_ng{ng}_ig{ig}{_ig_suffix}_sg{sg}"
nl_tmpl = "wrfout_d01_2018-07-05_120000_woinput_NL{}"
ad_dir  = f"./wrfout/woinput_absG_i0M_NOV0.02_sg_{ng}"
ad_tmpl = "wrfout_d01_2018-07-05_120000_woinput_AD{}"
num_members = 100
n_workers   = 6

# ==========================
# Load custom colormap
# ==========================
def load_rgb_cmap(filepath):
    colors = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 3:
                try:
                    r, g, b = int(parts[0]), int(parts[1]), int(parts[2])
                    colors.append((r / 255.0, g / 255.0, b / 255.0))
                except ValueError:
                    continue
    return LinearSegmentedColormap.from_list('testcmap', colors, N=len(colors))

custom_cmap = load_rgb_cmap("customcmap.rgb")


# ==========================
# Helper: parallel readers for Mode 2
# ==========================
def read_nl_rainnc(i):
    """Read terminal (last time step) RAINNC from NL member i."""
    fpath = os.path.join(nl_dir, nl_tmpl.format(i))
    try:
        with Dataset(fpath, 'r') as nc:
            return nc.variables['RAINNC'][-1, :, :]
    except Exception as e:
        print(f"[NL] Error member {i}: {e}")
        return None


def read_ad_rainnc(i):
    """Read t_index=6 RAINNC from AD member i."""
    fpath = os.path.join(ad_dir, ad_tmpl.format(i))
    try:
        with Dataset(fpath, 'r') as nc:
            return nc.variables['RAINNC'][6, :, :]
    except Exception as e:
        print(f"[AD] Error member {i}: {e}")
        return None


# ==========================
# Common plotting helper
# ==========================
def setup_map(ax, lons, lats):
    ax.set_extent([
        np.min(lons),
        np.max(lons) - 0.2,
        np.min(lats) + 0.1,
        np.max(lats) - 0.05,
    ], crs=ccrs.PlateCarree())
    ax.set_position([0.1, 0.1, 0.75, 0.75])
    ax.coastlines(resolution='10m', linewidth=1.0)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=0.3)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', alpha=0.8)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 18, 'color': 'black', 'weight': 'normal'}
    gl.ylabel_style = {'size': 18, 'color': 'black', 'weight': 'normal'}
    return gl


# ==========================
# 3x3 pseudo-Gaussian weighted mean
# ==========================
W = np.array([
    [0.25, 0.5, 0.25],
    [0.5,  1.0, 0.5 ],
    [0.25, 0.5, 0.25],
])

def pseudo_gaussian_weighted_mean_3x3(field, cy, cx, W):
    ny_f, nx_f = field.shape
    total_val = 0.0
    total_w   = 0.0
    for di in range(3):
        for dj in range(3):
            jj = cy + (di - 1)
            ii = cx + (dj - 1)
            if jj < 0 or jj >= ny_f or ii < 0 or ii >= nx_f:
                continue
            w = W[di, dj]
            total_val += w * field[jj, ii]
            total_w   += w
    return total_val / total_w if total_w > 0 else np.nan


# ==========================
# Load geographic info
# ==========================
with Dataset(geo_file_path, 'r') as nc:
    lats = nc.variables['XLAT'][0, :, :]
    lons = nc.variables['XLONG'][0, :, :]
    if mode == 1:
        QVAPOR_sfc = nc.variables['QVAPOR'][0, 0, :, :]
        nz = nc.variables['QVAPOR'].shape[1]

ny, nx = lats.shape
print(f"Grid: ny={ny}, nx={nx}")


# ================================================================
# Mode 1: QVAPOR change ratio
# ================================================================
if mode == 1:
    print("=== Mode 1: QVAPOR change ratio ===")

    flat_data = np.loadtxt(dat_path)
    expected_size = nz * ny * nx
    if flat_data.size != expected_size:
        raise ValueError(f"dat size {flat_data.size} != expected {expected_size}")

    A_QVAPOR_3D = flat_data.reshape((nz, ny, nx), order='F')
    A_QVAPOR_2D = A_QVAPOR_3D[0, :, :]

    # Normalize by Frobenius norm, then scale
    frob = np.linalg.norm(A_QVAPOR_2D, 'fro')
    print(f"Frobenius norm (raw): {frob:.6e}")
    A_QVAPOR_norm = A_QVAPOR_2D / frob * TARGET_NORM * ig
    print(f"Frobenius norm (after norm+scale by ig={ig}): {np.linalg.norm(A_QVAPOR_norm, 'fro'):.6e}")

    # Divide by surface QVAPOR; mask near-zero values
    mask = np.abs(QVAPOR_sfc) < 1e-10
    ratio = np.where(mask, np.nan, A_QVAPOR_norm / QVAPOR_sfc)
    print(f"Ratio range: [{np.nanmin(ratio):.6e}, {np.nanmax(ratio):.6e}]")

    # Version 1: original
    ratio_orig = ratio.copy()

    # Version 2: negatives only
    ratio_neg = ratio.copy()
    ratio_neg[ratio_neg > 0] = 0.0

    # Version 3: positives only
    ratio_pos = ratio.copy()
    ratio_pos[ratio_pos < 0] = 0.0

    # Version 4: top 8 most negative pixels (selected on A_QVAPOR_norm)
    ratio_top8 = np.zeros_like(ratio)
    neg_mask = A_QVAPOR_norm < 0
    neg_vals = A_QVAPOR_norm[neg_mask]
    if len(neg_vals) >= 8:
        threshold = np.sort(neg_vals)[7]  # 8th most negative
        keep_mask = neg_mask & (A_QVAPOR_norm <= threshold)
        ratio_top8[keep_mask] = ratio[keep_mask]
    else:
        ratio_top8[neg_mask] = ratio[neg_mask]


# ================================================================
# Mode 2: NL_mean - AD_mean RAINNC
# ================================================================
elif mode == 2:
    print("=== Mode 2: RAINNC diff (NL - AD) ===")

    members = list(range(1, num_members + 1))

    if ng == 0:
        # Single-member mode: ng=0 means no noise, only member 1 exists
        print("ng=0: single-member mode, reading member 1 only")
        nl_data = read_nl_rainnc(1)
        ad_data = read_ad_rainnc(1)
        if nl_data is None or ad_data is None:
            raise RuntimeError("Single-member read failed.")
        nl_mean = nl_data
        ad_mean = ad_data
    else:
        # ensemble average mode
        with Pool(processes=n_workers) as pool:
            nl_results = pool.map(read_nl_rainnc, members)
            ad_results = pool.map(read_ad_rainnc, members)

        nl_valid = [r for r in nl_results if r is not None]
        ad_valid = [r for r in ad_results if r is not None]
        print(f"NL valid: {len(nl_valid)}/{num_members}, AD valid: {len(ad_valid)}/{num_members}")

        if not nl_valid or not ad_valid:
            raise RuntimeError("Insufficient valid members for averaging.")

        nl_mean = np.mean(np.stack(nl_valid, axis=0), axis=0)
        ad_mean = np.mean(np.stack(ad_valid, axis=0), axis=0)

    diff = nl_mean - ad_mean
    print(f"RAINNC diff range: [{diff.min():.4f}, {diff.max():.4f}]")

    # ---- plot params ----
    field     = diff
    abs_max   = np.max(np.abs(field))
    vmin, vmax = -abs_max, abs_max
    cb_label  = 'RAINNC difference (mm)'
    title_str = 'RAINNC: NL(terminal) − AD(t=6) ensemble average'
    save_name = f"./figs/RAINNC_diff_NL_AD_ng{ng}_ig{ig}{_ig_suffix}_sg{sg}.png"

else:
    raise ValueError(f"Unknown mode: {mode}. Use 1 or 2.")


# ==========================
# Plot
# ==========================
mpl.rcParams['font.family'] = 'Helvetica'

if mode == 1:
    plots = [
        (ratio_orig, f'QVAPOR change ratio (sg={sg})',                       f'ratio_AQVAPOR_QVAPOR_sg{sg}'),
        (ratio_neg,  f'QVAPOR change ratio — negatives only (sg={sg})',       f'ratio_AQVAPOR_QVAPOR_sg{sg}_negonly'),
        (ratio_pos,  f'QVAPOR change ratio — positives only (sg={sg})',       f'ratio_AQVAPOR_QVAPOR_sg{sg}_posonly'),
        (ratio_top8, f'QVAPOR change ratio — top 8 negative pixels (sg={sg})',f'ratio_AQVAPOR_QVAPOR_sg{sg}_top8neg'),
    ]
    for data, title_str, fname in plots:
        fig = plt.figure(figsize=(10, 8))
        ax  = plt.axes(projection=ccrs.PlateCarree())
        setup_map(ax, lons, lats)

        mesh = ax.pcolormesh(lons, lats, data, cmap='seismic',
                             vmin=-1, vmax=1, shading='auto')
        cb = plt.colorbar(mesh, orientation='vertical', shrink=0.85, pad=0.04, aspect=30, ax=ax)
        cb.set_label('QVAPOR change ratio', fontsize=16)
        cb.ax.tick_params(labelsize=18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.set_title(title_str, fontweight='bold', pad=10, fontsize=20)
        draw_small_box(ax, lons, lats, x_center=16, y_center=34, width=4)

        save_name = f"./figs/{fname}.png"
        plt.savefig(save_name, dpi=600, bbox_inches='tight')
        plt.close()
        print(f"Saved: {save_name}")

else:
    fig = plt.figure(figsize=(10, 8))
    ax  = plt.axes(projection=ccrs.PlateCarree())
    setup_map(ax, lons, lats)

    mesh = ax.pcolormesh(lons, lats, field, cmap=custom_cmap,
                         vmin=-20, vmax=20, shading='auto')
    cb = plt.colorbar(mesh, orientation='vertical', shrink=0.85, pad=0.04, aspect=30, ax=ax)
    cb.set_label(cb_label, fontsize=16)
    cb.ax.tick_params(labelsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.set_title(title_str, fontweight='bold', pad=10, fontsize=20)
    draw_small_box(ax, lons, lats, x_center=16, y_center=34, width=4)

    x_center, y_center = 16, 34
    cy, cx = y_center, x_center
    print(f"Box center: cy={cy}, cx={cx}, lon={lons[cy,cx]:.3f}, lat={lats[cy,cx]:.3f}")
    wm = pseudo_gaussian_weighted_mean_3x3(field, cy, cx, W)
    print(f"3x3 pseudo-Gaussian weighted mean at box center: {wm:.4f}")

    plt.savefig(save_name, dpi=600, bbox_inches='tight')
    plt.show()
    print(f"Figure saved to {save_name}")