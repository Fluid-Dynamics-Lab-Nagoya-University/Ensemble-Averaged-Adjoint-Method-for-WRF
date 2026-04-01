#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots color contour for precipitation reduction
of different parameter combinations.

by Shan Jiang, FDL, Nagoya University

"""

import os
import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

from func_analysis import compute_target_rr_custom_kernel
from plotting_tools import plot_contour

# ================= Parameters =================
base_path = "./wrfout"
t_idx = 6
center_x_Gker = 34 #34 for absG_i0M, 35 for x30y36
center_y_Gker = 16 #16 for absG_i0M, 29 for x30y36
nCases = 100

ig_list   = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]
#ig_list   = ["0.1", "0.2", "0.3", "0.4"]
#ng_values = ["0.1"]
ng_values = ["0"]
#ng_values = ["0.1", "0.2", "0.3", "0.4", "0.5"]
s_values = ["0", "0.001", "0.01", "0.05", "0.1", "0.15", "0.2", "0.3", "0.4", "0.5"]  # sg values; equals ng when use_ng_eq_sg=True
#s_values = ["0"]

MAX_WORKERS = 8  # adjust based on disk I/O and CPU availability
input_mode = 1  # 1: no suffix (both signs); 2: positive input only (P); 3: negative input only (N)
actr_mode = 0   # 1: include "actr" in directory name; 0: omit

s_mode = 1  # 1: sg 2:sv
n_mode = 1  # 1: ng 2: nv

run_mode = 2   # 1: read NetCDF, compute, save npy, and plot; 2: load existing npy and plot only
eval_mode = 2  # 1: single-point evaluation at center; 2: 3x3 pseudo-Gaussian weighted average

# ======================================================
use_ng_eq_sg = False   # True: ng equals sg for each iteration; False: ng fixed at "0"

# ======================================================
s_tag = "sg" if s_mode == 1 else "sv"
n_tag = "ng" if n_mode == 1 else "nv"

ref_s_tag = "sg" if n_mode == 1 else "sv" # because in RAINNCdiff calculation, AD case is also used.
actr_tag = "actr" if actr_mode == 1 else ""

def make_save_tag(ng):
    suffix = "" if input_mode == 1 else ("P" if input_mode == 2 else "N")
    actr_label = "actr" if actr_mode == 1 else ""
    s_label = "sg" if s_mode == 1 else "sv"
    eval_tag = "" if eval_mode == 1 else "_G3x3"

    if use_ng_eq_sg:
        ng_part = "ngEQsg"
    else:
        ng_part = f"{n_tag}{ng}"

    return f"absG_i0M_NOV0.02_{ng_part}{suffix}{actr_label}_{s_label}{eval_tag}"


# ================= Top-level functions for subprocess pickling =================
def _read_rainnc_at_t(path, t_idx):
    try:
        with Dataset(path, 'r') as ds:
            return ds.variables["RAINNC"][t_idx, :, :].astype(np.float32, copy=False)
    except FileNotFoundError:
        return None
    except Exception as e:
        print(f"[ERROR] {path}: {e}")
        return None

def mean_over_cases(paths, desc=None, max_workers=MAX_WORKERS):
    total = None
    count = 0
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_read_rainnc_at_t, p, t_idx) for p in paths]
        it = futures
        if desc:
            it = tqdm(futures, total=len(futures), desc=desc, leave=False)
        for fut in it:
            arr = fut.result()
            if arr is None:
                continue
            if total is None:
                total = arr.copy()
            else:
                total += arr
            count += 1
    if count == 0:
        raise RuntimeError("No valid cases found for averaging.")
    return total / np.float32(count)

def build_ref_paths(ng):
    folder = f"{base_path}/woinput_absG_i0M_NOV0.02_{ref_s_tag}_{ng}"
    num_cases = 1 if ng == "0" else nCases
    return [f"{folder}/wrfout_d01_2018-07-05_120000_woinput_AD{i}" for i in range(1, num_cases + 1)]

def build_ad_paths(ng, ig, sg, input_mode):
    suffix = "" if input_mode == 1 else ("P" if input_mode == 2 else "N")
    folder = f"{base_path}/absG_i0M_NOV0.02_{n_tag}{ng}_ig{ig}{suffix}{actr_tag}_{s_tag}{sg}"
    num_cases = 1 if ng == "0" else nCases
    return [f"{folder}/wrfout_d01_2018-07-05_120000_woinput_NL{i}" for i in range(1, num_cases + 1)]


# ================= Main =================
def main():
    os.makedirs("./dats", exist_ok=True)
    os.makedirs("./figs", exist_ok=True)

    save_tag = make_save_tag(ng_values[0])

    # ---------- run_mode == 2: load existing npy and plot ----------
    if run_mode == 2:
        deltaP_mat = np.load(f"./dats/deltaP_mat_{save_tag}.npy")
        print(f"[INFO] Loaded npy files for tag: {save_tag}")

    # ---------- run_mode == 1: read NetCDF and compute ----------
    else:
        deltaP_mat = np.zeros((len(s_values), len(ig_list)), dtype=np.float32)

        # Precompute reference fields.
        # When use_ng_eq_sg=True, ng varies with sg, so precompute a ref for each sg value.
        # When use_ng_eq_sg=False, ng is fixed at "0", so only one ref is needed.
        ref_keys = s_values if use_ng_eq_sg else ng_values
        ref_cache = {}
        for ng_key in tqdm(ref_keys, desc="Precompute ref per ng"):
            ref_paths = build_ref_paths(ng_key)
            ref_cache[ng_key] = mean_over_cases(ref_paths, desc=f"  ref ng={ng_key}", max_workers=MAX_WORKERS)

        for j_s, s in enumerate(tqdm(s_values, desc=f"outer layer {s_tag}")):
            ng = s if use_ng_eq_sg else ng_values[0]
            ref_mean_rainnc = ref_cache[ng]
            for i_ig, ig in enumerate(tqdm(ig_list, desc=f"  -- inner ig@{s_tag}={s}", leave=False)):
                ad_paths = build_ad_paths(ng, ig, s, input_mode)
                mean_rainnc = mean_over_cases(ad_paths, desc="    AD cases", max_workers=MAX_WORKERS)

                if eval_mode == 1:
                    _, _, _, targetP = compute_target_rr_custom_kernel(
                        ref_mean_rainnc, mean_rainnc, center_x_Gker, center_y_Gker
                    )
                else:
                    # 3x3 pseudo-Gaussian weighted average centered at (center_x_Gker, center_y_Gker)
                    weights = np.array([
                        [0.25, 0.5, 0.25],
                        [0.5,  1.0, 0.5 ],
                        [0.25, 0.5, 0.25]
                    ], dtype=np.float32)

                    weighted_sum = 0.0
                    weight_sum   = 0.0

                    for j_dy, dy in enumerate([-1, 0, 1]):
                        for j_dx, dx in enumerate([-1, 0, 1]):
                            w = float(weights[j_dy, j_dx])
                            if w == 0.0:
                                continue

                            cx = center_x_Gker + dx
                            cy = center_y_Gker + dy

                            _, _, _, tp = compute_target_rr_custom_kernel(
                                ref_mean_rainnc, mean_rainnc, cx, cy
                            )
                            weighted_sum += w * tp
                            weight_sum   += w

                    targetP = weighted_sum / weight_sum

                deltaP_mat[j_s, i_ig] = targetP

                suffix = "" if input_mode == 1 else ("P" if input_mode == 2 else "N")
                print(f"[{n_tag}={ng} | ig={ig}{suffix}{actr_tag} | {s_tag}={s}]  ΔP = {targetP:.3f} mm")

        np.save(f"./dats/deltaP_mat_{save_tag}.npy", deltaP_mat)
        print(f"[INFO] Saved: ./dats/deltaP_mat_{save_tag}.npy")

    # ---------- Plot ----------
    plot_contour(
        ig_list, s_values, deltaP_mat,
        title=f"$-\Delta P$ (mm) vs ig and {s_tag}",
        label="$-\Delta P$ (mm)",
        filename=f"./figs/deltaP_contour_{save_tag}.png",
        vmin=0, vmax=10.5
    )

    return deltaP_mat  # for workspace


# ================= Entry point =================
if __name__ == "__main__":
    deltaP_mat = main()
