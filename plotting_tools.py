#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A helper function for plotting.

by Shan Jiang, FDL, Nagoya University
"""

import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt


def draw_small_box(ax, lons, lats, x_center=30, y_center=36, width=4, color='magenta'):
    """
    Draw a rectangle on a map centered at (x_center, y_center) in grid index space.

    Parameters:
        ax: matplotlib/cartopy axis object
        lons, lats: 2D longitude and latitude arrays
        x_center, y_center: grid center index
        width: line width
        color: edge color
    """
    x0 = x_center - 2
    x1 = x_center + 2
    y0 = y_center - 2
    y1 = y_center + 2

    lon_min, lat_min = lons[y0, x0], lats[y0, x0]
    lon_max, lat_max = lons[y1, x1], lats[y1, x1]

    rect = mpatches.Rectangle(
        (lon_min, lat_min),
        lon_max - lon_min,
        lat_max - lat_min,
        linewidth=width,
        edgecolor=color,
        facecolor='none',
        transform=ccrs.PlateCarree(),
        zorder=10
    )
    ax.add_patch(rect)


def plot_contour(iv_list, sv_list, data_matrix, title, label, filename,
                 vmin=None, vmax=None):
    """
    Plot a discrete heatmap (no interpolation) using imshow.
    x-axis: iv, y-axis: sv, color: data_matrix values.

    vmin, vmax: color scale range; defaults to data_matrix min/max if not specified.
    """

    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['font.size'] = 18

    if vmin is None:
        vmin = float(np.nanmin(data_matrix))
    if vmax is None:
        vmax = float(np.nanmax(data_matrix))

    fig, ax = plt.subplots(figsize=(8, 6))

    im = ax.imshow(
        data_matrix, cmap='jet', origin='lower', aspect='auto',
        vmin=vmin, vmax=vmax
    )

    ax.set_xticks(np.arange(len(iv_list)))
    ax.set_yticks(np.arange(len(sv_list)))
    ax.set_xticklabels(iv_list, fontsize=18)
    ax.set_yticklabels(sv_list, fontsize=18)

    ax.set_xlabel("Input Rate", fontsize=20)
    ax.set_ylabel("Noise for calculating sensitivity")
    ax.set_title(title, fontsize=20, fontweight='bold')

    # Value labels
    for i in range(len(sv_list)):
        for j in range(len(iv_list)):
            value = data_matrix[i, j]
            ax.text(j, i, f"{value:.2f}",
                    ha='center', va='center', color='white', fontsize=14)

    cbar = plt.colorbar(im)
    cbar.set_label(label, fontsize=18)
    cbar.ax.tick_params(labelsize=18)

    plt.tight_layout()
    plt.savefig(filename, dpi=600)
    plt.close()
