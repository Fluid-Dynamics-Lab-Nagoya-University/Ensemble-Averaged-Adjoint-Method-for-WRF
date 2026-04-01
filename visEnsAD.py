#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ensemble averaging and visualizing AD-calculated sensitivity.

by Shan Jiang, FDL, Nagoya University
"""

import numpy as np
from netCDF4 import Dataset
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from plotting_tools import draw_small_box

# ==========================
# Configuration
# ==========================

param_name  = "sg"   # change to "sv" to switch
param_value = 0.3

base_dir = f"./wrfout/woinput_absG_i0M_NOV0.02_{param_name}_{param_value}"
file_prefix = "wrfout_d01_2018-07-05_120000_woinput_AD"
num_members = 100
t_index = -1

# ---- control options ----
normalize_flag = 1   # 0: no normalization, 1: normalize to [-1, 1]
neg_top_flag   = 1   # 0: normal plot, 1: normalize then keep only the 8 largest-magnitude negative points

# ==========================
# Load ensemble data
# ==========================
A_QVAPOR_list = []
for i in range(1, num_members + 1):
    filename = os.path.join(base_dir, f"{file_prefix}{i}")
    if not os.path.isfile(filename):
        print(f"Warning: File not found: {filename} (skipped)")
        continue
    try:
        with Dataset(filename, 'r') as nc:
            qvapor = nc.variables['A_QVAPOR'][t_index, :, :, :]
            A_QVAPOR_list.append(qvapor)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        continue

if not A_QVAPOR_list:
    raise RuntimeError("No A_QVAPOR data was successfully read.")

# ==========================
# Compute ensemble average
# ==========================
A_QVAPOR_stack = np.stack(A_QVAPOR_list, axis=0)
A_QVAPOR_mean = np.mean(A_QVAPOR_stack, axis=0)

os.makedirs("./dats", exist_ok=True)
os.makedirs("./figs", exist_ok=True)

np.savetxt(f"./dats/A_QVAPOR_mean_absG_i0M_NOV0.02_{param_name}{param_value}.dat",
           A_QVAPOR_mean.reshape(-1, order='F'), fmt="%.14e")
print("A_QVAPOR ensemble average saved.")

# ==========================
# Load geographic info
# ==========================
with Dataset('./clean/wrfinput_d01_clean', 'r') as nc:
    lats = nc.variables['XLAT'][0, :, :]
    lons = nc.variables['XLONG'][0, :, :]

layer_index = 0
A_QVAPOR_layer = A_QVAPOR_mean[layer_index, :, :].copy()

if A_QVAPOR_layer.shape != lats.shape:
    raise ValueError(f"Grid shape mismatch: A_QVAPOR_mean shape = {A_QVAPOR_mean.shape}, "
                     f"but lat/lon shape = {lats.shape}")

# ==========================
# Data processing
# ==========================
if neg_top_flag == 1:
    # Step 1: normalize
    abs_max_all = np.max(np.abs(A_QVAPOR_layer))
    if abs_max_all == 0:
        raise ValueError("Max absolute value is zero; cannot normalize before neg_top selection.")
    norm_layer = A_QVAPOR_layer / abs_max_all

    # Step 2: keep only 8 strongest negative points
    flat = norm_layer.flatten(order='F')
    neg_idx = np.where(flat < 0)[0]
    if neg_idx.size > 0:
        k = min(8, neg_idx.size)
        pick_rel = np.argsort(np.abs(flat[neg_idx]))[-k:]
        keep_idx = neg_idx[pick_rel]
        kept = np.zeros_like(flat)
        kept[keep_idx] = flat[keep_idx]
        A_QVAPOR_layer = kept.reshape(norm_layer.shape, order='F')
        print(f"neg_top_flag applied: kept {k} negative points after normalization.")
    else:
        A_QVAPOR_layer = np.zeros_like(norm_layer)
        print("Warning: no negative values found after normalization; map is all zeros.")

    vmin, vmax = -1.0, 1.0
    suffix = "_negTop_norm"

else:
    if normalize_flag == 1:
        abs_max = np.max(np.abs(A_QVAPOR_layer))
        if abs_max == 0:
            raise ValueError("Max absolute value is zero; cannot normalize.")
        A_QVAPOR_layer = A_QVAPOR_layer / abs_max
        vmin, vmax = -1.0, 1.0
        suffix = "_norm"
    else:
        vmin, vmax = -6500, 6500
        suffix = ""

# ==========================
# Plotting
# ==========================
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())

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

mesh = ax.pcolormesh(lons, lats, A_QVAPOR_layer, cmap='seismic',
                     vmin=vmin, vmax=vmax, shading='auto')

mpl.rcParams['font.family'] = 'Helvetica'
cb = plt.colorbar(mesh, orientation='vertical', shrink=0.85, pad=0.04, aspect=30, ax=ax)
cb.set_label('A_QVAPOR', fontsize=20)
cb.ax.tick_params(labelsize=18)
ax.tick_params(axis='both', which='major', labelsize=18)

ax.set_title(
    f'Sensitivity calculated with {param_value} noise'
    + (" (normalized; top 8 negative)" if neg_top_flag == 1 else
       (" (normalized)" if normalize_flag == 1 else "")),
    fontweight='bold', pad=10, fontsize=20
)

draw_small_box(ax, lons, lats, x_center=16, y_center=34, width=4)

plt.savefig(f"./figs/A_QVAPORabsG_i0M_NOV0.02_ens_{param_name}{param_value}{suffix}.png",
            dpi=600, bbox_inches='tight')
plt.show()
