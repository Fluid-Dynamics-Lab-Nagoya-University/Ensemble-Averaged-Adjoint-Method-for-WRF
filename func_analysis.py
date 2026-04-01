#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A helper function for analyzing precipitation difference
at designated Gaussian area.

by Shan Jiang, FDL, Nagoya University

"""

import numpy as np

def compute_target_rr_custom_kernel(rain_ref, rain_ctrl, center_x, center_y):
    """
    Compute rainfall reduction at the target grid point using a single-point kernel.
    The 3x3 pseudo-Gaussian weighted evaluation is handled externally by the caller.

    Parameters:
    - rain_ref:   2D numpy array, reference simulation rainfall
    - rain_ctrl:  2D numpy array, control simulation rainfall
    - center_x, center_y: center grid point index (int)

    Returns:
    - P_ref_weighted:  weighted rainfall in reference case
    - P_ctrl_weighted: weighted rainfall in control case
    - RR: rainfall reduction rate (nan if P_ref_weighted == 0)
    - RV: rainfall reduction value (nan if P_ref_weighted == 0)
    """

    assert rain_ref.shape == rain_ctrl.shape, "Input arrays must have the same shape!"

    # Single-point kernel: reads only the center value.
    # Alternative 3x3 pseudo-Gaussian kernel (not used here; weighting done externally):
    #   [0.025, 0.05,  0.025]
    #   [0.05,  0.1,   0.05 ]
    #   [0.025, 0.05,  0.025]
    kernel = np.array([
        [0, 0, 0],
        [0, 1, 0],
        [0, 0, 0]
    ])

    # Define 3x3 window boundaries
    x_start = center_x - 1
    x_end   = center_x + 2
    y_start = center_y - 1
    y_end   = center_y + 2

    # Check bounds
    if x_start < 0 or y_start < 0 or x_end > rain_ref.shape[0] or y_end > rain_ref.shape[1]:
        raise ValueError("Target region is too close to the boundary for a 3x3 kernel.")

    # Extract patches
    ref_patch = rain_ref[x_start:x_end, y_start:y_end]
    ctrl_patch = rain_ctrl[x_start:x_end, y_start:y_end]

    # Compute weighted rainfall
    P_ref_weighted = np.sum(ref_patch * kernel)
    P_ctrl_weighted = np.sum(ctrl_patch * kernel)

    # Compute reduction rate
    if P_ref_weighted == 0:
        RR = np.nan
        RV = np.nan
    else:
        RR = (P_ref_weighted - P_ctrl_weighted) / P_ref_weighted
        RV = P_ref_weighted - P_ctrl_weighted

    return P_ref_weighted, P_ctrl_weighted, RR, RV

