#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metric phi = (1 - cos(mean_s, s_det)) * gamma
  - deflection term: how much ensemble mean deviates from deterministic gradient
  - Coherence term: gamma = ||mean_s|| / mean(||s_k||)

Also tries variants:
  phi2 = (1 - cos) * gamma^2
  phi4 = sin(angle) * ||mean_s||

Input: dataset-root wrfout/diagnostics/, overridden by WRF_DIAGNOSTIC_ROOT.
Switch modes (applied AFTER ensemble averaging):
  1 = unfiltered
  2 = keep only negative values (positive -> 0)
  3 = keep only positive values (negative -> 0)
  4 = keep top 8 most negative points only

NOTE: coherence (gamma) is always computed from the FULL (unfiltered) gradient
fields, since it characterizes the statistical consistency of the ensemble
itself, independent of actuation strategy. The switch filtering is applied
only to the ensemble-averaged sensitivity and the deterministic gradient
when computing deflection and the composite metrics.

8-core parallel + tqdm.

by Shan Jiang, FDL, Nagoya University
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Pool
from tqdm import tqdm
import os
from pathlib import Path

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Helvetica']
mpl.rcParams['font.size'] = 14

# Parameters
base_root = os.environ.get('WRF_DIAGNOSTIC_ROOT', str(Path(__file__).resolve().parents[1] / 'wrfout' / 'diagnostics'))
sg_list = [0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
n_members = 100
n_workers = 8

# ===== SWITCH =====
# 1: unfiltered
# 2: keep only negative values (positive -> 0)
# 3: keep only positive values (negative -> 0)
# 4: keep top 8 most negative points, rest -> 0
switch = 1

# Font size offset: adjusts ALL font sizes by this amount (in pt)
font_size_offset = 1

# Unified figure and plot area size
fig_width, fig_height = 9, 5.5
plot_rect = [0.12, 0.14, 0.82, 0.76]

# Fixed accumulated-precipitation-reduction references for the plot overlay.
best_ctrl = {
    0:     1.8,
    0.001: 1.32,
    0.01:  1.60,
    0.05:  4.43,
    0.1:   5.77,
    0.15:  8.24,
    0.2:   8.86,
    0.3:   10.17,
    0.4:   9.02,
    0.5:   6.96,
}

# Output directory
fig_dir = "./figs"
os.makedirs(fig_dir, exist_ok=True)

# Switch label for filenames and titles
switch_labels = {
    1: "all",
    2: "neg_only",
    3: "pos_only",
    4: "top8neg",
}
switch_suffix = switch_labels.get(switch, f"sw{switch}")

# Font size helper
def fs(base_size):
    """Apply font_size_offset to a base font size."""
    return base_size + font_size_offset

# Field filtering (applied AFTER ensemble averaging)
def apply_switch(field_2d, switch_mode):
    """Apply switch-based filtering to a 2D field.

    In the corrected version, this is applied to the ensemble-averaged
    sensitivity (or deterministic gradient), NOT to individual members.
    """
    out = field_2d.copy()
    if switch_mode == 1:
        pass  # no filtering
    elif switch_mode == 2:
        out[out > 0] = 0.0
    elif switch_mode == 3:
        out[out < 0] = 0.0
    elif switch_mode == 4:
        out[out > 0] = 0.0
        flat = out.ravel()
        neg_indices = np.where(flat < 0)[0]
        if len(neg_indices) > 8:
            sorted_neg = neg_indices[np.argsort(flat[neg_indices])]
            keep = sorted_neg[:8]
            mask = np.zeros_like(flat)
            mask[keep] = flat[keep]
            out = mask.reshape(out.shape)
    else:
        raise ValueError(f"Unknown switch value: {switch_mode}")
    return out

# Helpers
def get_filepath(sg, member_id):
    if sg == 0:
        sg_str = "0"
    else:
        sg_str = f"{sg:g}"
    folder = f"woinput_absG_i0M_NOV0.02_sg_{sg_str}"
    fname = f"wrfout_d01_2018-07-05_120000_woinput_AD{member_id}"
    return f"{base_root}/{folder}/{fname}"

def _read_one_raw(args):
    """Read raw (unfiltered) gradient field for one member."""
    sg, member_id = args
    fpath = get_filepath(sg, member_id)
    try:
        with Dataset(fpath, 'r') as nc:
            grad_2d = nc.variables['A_QVAPOR'][-1, 0, :, :]
        return grad_2d.astype(np.float64)
    except Exception as e:
        raise RuntimeError(f"Failed to read gradient member: {fpath}") from e

def read_all_members_raw(sg, pool):
    """Read all members' raw (unfiltered) gradient fields."""
    n_mem = 1 if sg == 0 else n_members
    tasks = [(sg, k) for k in range(1, n_mem + 1)]
    results = list(tqdm(
        pool.imap(_read_one_raw, tasks),
        total=n_mem,
        desc=f"  sg={sg:<6g}",
        ncols=80,
        leave=True
    ))
    return results

# Main
if __name__ == '__main__':
    print(f"Switch mode: {switch} ({switch_suffix})")
    print("NOTE: filtering is applied AFTER ensemble averaging.")
    print("      coherence (gamma) is computed from FULL (unfiltered) fields.\n")

    # Step 1: Read deterministic gradient (sg=0), FULL field
    print("Reading deterministic gradient (sg=0)...")
    fpath_det = get_filepath(0, 1)
    try:
        with Dataset(fpath_det, 'r') as nc:
            g_det_2d_full = nc.variables['A_QVAPOR'][-1, 0, :, :].astype(np.float64)
        print(f"  ||g_det (full)|| = {np.linalg.norm(g_det_2d_full):.6e}")
    except Exception as e:
        raise RuntimeError(f"Failed to read deterministic gradient: {fpath_det}") from e

    # Apply switch to deterministic gradient for deflection computation
    g_det_2d_filtered = apply_switch(g_det_2d_full, switch)
    g_det_filtered = g_det_2d_filtered.ravel()
    norm_det_filtered = np.linalg.norm(g_det_filtered)
    print(f"  ||g_det (filtered, {switch_suffix})|| = {norm_det_filtered:.6e}")

    # Step 2: For each sg, read ALL members raw, then:
    #   - compute coherence from FULL fields
    #   - compute ensemble mean from FULL fields
    #   - apply switch to ensemble mean
    #   - compute deflection from filtered fields
    results = {}

    with Pool(processes=n_workers) as pool:
        for sg in sg_list:
            if sg == 0:
                continue
            print(f"\n--- C_pert = {sg} ---")
            grads_2d = read_all_members_raw(sg, pool)

            # --- Coherence: computed from FULL (unfiltered) fields ---
            grads_full_flat = np.array([g.ravel() for g in grads_2d])
            norms_full = np.linalg.norm(grads_full_flat, axis=1)
            norms_full_safe = np.clip(norms_full, 1e-30, None)

            mean_grad_full = np.mean(grads_full_flat, axis=0)
            norm_mean_full = np.linalg.norm(mean_grad_full)
            mean_norm_full = np.mean(norms_full_safe)

            gamma = norm_mean_full / mean_norm_full

            # --- Deflection: computed from FILTERED fields ---
            # Ensemble mean (2D) -> apply switch -> flatten
            mean_grad_2d = np.mean(np.array(grads_2d), axis=0)
            mean_grad_2d_filtered = apply_switch(mean_grad_2d, switch)
            mean_grad_filtered = mean_grad_2d_filtered.ravel()
            norm_mean_filtered = np.linalg.norm(mean_grad_filtered)

            if norm_mean_filtered > 1e-30 and norm_det_filtered > 1e-30:
                cos_det = np.dot(mean_grad_filtered, g_det_filtered) / (norm_mean_filtered * norm_det_filtered)
                cos_det = np.clip(cos_det, -1.0, 1.0)
            else:
                cos_det = 1.0  # no signal

            deflection = 1.0 - cos_det
            sin_det = np.sqrt(1.0 - cos_det**2)

            # --- Composite metrics ---
            phi1 = deflection * gamma
            phi2 = deflection * gamma**2
            phi4 = sin_det * norm_mean_filtered

            results[sg] = {
                'gamma': gamma,
                'cos_det': cos_det,
                'deflection': deflection,
                'norm_mean_full': norm_mean_full,
                'norm_mean_filtered': norm_mean_filtered,
                'phi1': phi1,
                'phi2': phi2,
                'phi4': phi4,
            }
            print(f"  gamma={gamma:.4f} (from full field)")
            print(f"  cos_det={cos_det:.4f} (from filtered field, {switch_suffix})")
            print(f"  phi1={phi1:.6f}, phi2={phi2:.6f}, phi4={phi4:.6e}")

    # Prepare arrays
    sg_ens = sorted(results.keys())
    sg_arr = np.array(sg_ens)
    gamma_arr = np.array([results[s]['gamma'] for s in sg_ens])
    deflection_arr = np.array([results[s]['deflection'] for s in sg_ens])
    phi1_arr = np.array([results[s]['phi1'] for s in sg_ens])
    phi2_arr = np.array([results[s]['phi2'] for s in sg_ens])
    phi4_arr = np.array([results[s]['phi4'] for s in sg_ens])

    sg_all = np.array(sg_list)
    ctrl_arr = np.array([best_ctrl.get(s, np.nan) for s in sg_list])

    # PLOT 1: deflection and Gamma separately
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax1 = fig.add_axes(plot_rect)

    ax1.plot(sg_arr, deflection_arr, 'o-', color='#FFCC22', lw=2, ms=7,
             label=r'$1 - \cos(\overline{\mathbf{s}},\, \mathbf{s})$ (deflection)')
    ax1.plot(sg_arr, gamma_arr, 's-', color='#2166AC', lw=2, ms=7,
             label=r'$\gamma$ (coherence)')
    ax1.set_xlabel(r'$C_{\mathrm{pert}}$', fontsize=fs(16))
    ax1.set_ylabel('Value', fontsize=fs(14))
    ax1.axvline(x=0.3, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax1.legend(fontsize=fs(12))
    ax1.set_title(f'Deflection vs Coherence: competing trends [{switch_suffix}]',
                  fontweight='bold', pad=12, fontsize=fs(14))
    ax1.tick_params(labelsize=fs(12))
    fig.savefig(f"{fig_dir}/plot1_deflection_coherence_{switch_suffix}.png",
                dpi=600, bbox_inches='tight')
    plt.show()

    # PLOT 2: metric variants vs control performance (all normalized)
    def norm01(a):
        return (a - a.min()) / (a.max() - a.min() + 1e-30)

    ctrl_norm = (ctrl_arr - ctrl_arr.min()) / (ctrl_arr.max() - ctrl_arr.min() + 1e-30)

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_axes(plot_rect)

    ax.plot(sg_arr, norm01(phi1_arr), 'o-', lw=2, ms=6,
            label=r'$[1 - \cos(\overline{\mathbf{s}},\, \mathbf{s})]\,\gamma$')
    ax.plot(sg_arr, norm01(phi2_arr), 's-', lw=2, ms=6,
            label=r'$[1 - \cos(\overline{\mathbf{s}},\, \mathbf{s})]\,\gamma^2$')
    ax.plot(sg_all, ctrl_norm, 'k*--', lw=2.5, ms=12, alpha=0.8,
            label=r'$-\Delta\mathbf{P}^f$ (control performance)')

    ax.axvline(x=0.3, color='gray', ls=':', lw=1.5, alpha=0.7)
    ax.set_xlabel(r'$C_{\mathrm{pert}}$', fontsize=fs(16))
    ax.set_ylabel('Normalized value [0, 1]', fontsize=fs(14))
    ax.legend(fontsize=fs(12), ncol=2)
    ax.set_title(f'Metrics (normalized) vs control performance [{switch_suffix}]',
                 fontweight='bold', pad=12, fontsize=fs(14))
    ax.tick_params(labelsize=fs(12))
    fig.savefig(f"{fig_dir}/plot2_phi_normalized_{switch_suffix}.png",
                dpi=600, bbox_inches='tight')
    plt.show()

    # Print summary table
    print("\n" + "="*90)
    print(f"Switch mode: {switch} ({switch_suffix})")
    print(f"  - coherence (gamma): from FULL field")
    print(f"  - deflection & phi:  from FILTERED field ({switch_suffix})")
    print(f"{'C_pert':>8} {'gamma':>8} {'cos_det':>8} {'defl':>8} "
          f"{'phi1':>10} {'phi2':>10} {'phi4':>12} {'ctrl':>8}")
    print("-"*90)
    for sg in sg_ens:
        r = results[sg]
        c = best_ctrl.get(sg, np.nan)
        print(f"{sg:>8g} {r['gamma']:>8.4f} {r['cos_det']:>8.4f} {r['deflection']:>8.4f} "
              f"{r['phi1']:>10.6f} {r['phi2']:>10.6f} {r['phi4']:>12.6e} {c:>8.1f}")
    print("="*90)

    print("\nDone.")
