#!/usr/bin/env python3
"""Quantify finite-ensemble sampling uncertainty from existing cached WRF output.

This cache-only workflow never opens or runs WRF and never writes to the source cache.
It analyzes raw member-wise gradients and paired baseline/controlled objectives using
50-member sampling without replacement, subset convergence, and leave-one-out diagnostics.

by Shan Jiang, FDL, Nagoya University
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_CACHE_ROOT = HERE.parent / "adjoint_prediction_validation" / "cache"
DEFAULT_OUTPUT_DIR = HERE / "results" / "subsample_50"
DEFAULT_ANALYSIS_CACHE_ROOT = HERE / "cache_subsample_50"
SAMPLE_SIZE = 50
REQUIRED_CACHE_FILES = (
    "s_raw.npy",
    "j_baseline_mm.npy",
    "j_controlled_mm.npy",
    "metadata.json",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use the existing 100-member cache to quantify finite-ensemble sampling "
            "uncertainty; no WRF or NetCDF files are opened."
        )
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Complete adjoint-prediction cache; default: auto-discover one matching cache.",
    )
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--subsample-replicates", type=int, default=5000)
    parser.add_argument("--subset-repeats", type=int, default=1000)
    parser.add_argument(
        "--subset-sizes",
        type=int,
        nargs="+",
        default=[5, 10, 20, 30, 50, 75, 100],
    )
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--highlight-cpert", type=float, default=0.3)
    parser.add_argument(
        "--selected-cact",
        type=float,
        nargs="+",
        default=[0.3, 0.6, 0.9],
        help="Cact panels shown in the response percentile-interval figure.",
    )
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--analysis-cache-root", type=Path, default=DEFAULT_ANALYSIS_CACHE_ROOT)
    parser.add_argument("--from-cache", action="store_true", help="Require a matching analysis cache.")
    parser.add_argument("--overwrite", action="store_true", help="Replace results in the selected subsampling output directory.")
    parser.add_argument("--self-check", action="store_true", help="Run a small synthetic sampling check.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)

    if args.subsample_replicates < 100:
        parser.error("--subsample-replicates must be at least 100")
    if args.subset_repeats < 100:
        parser.error("--subset-repeats must be at least 100")
    if args.batch_size <= 0:
        parser.error("--batch-size must be positive")
    if any(size <= 0 for size in args.subset_sizes):
        parser.error("all --subset-sizes must be positive")
    return args


def cache_is_complete(path: Path) -> bool:
    if not path.is_dir():
        return False
    if not all((path / name).is_file() for name in REQUIRED_CACHE_FILES):
        return False
    try:
        metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return bool(metadata.get("complete"))


def resolve_cache(args: argparse.Namespace) -> Path:
    if args.cache_dir is not None:
        path = args.cache_dir.expanduser().resolve()
        if not cache_is_complete(path):
            raise FileNotFoundError(f"Cache is missing or incomplete: {path}")
        return path

    root = args.cache_root.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Cache root does not exist: {root}")
    candidates = [p for p in root.glob("cache_*") if cache_is_complete(p)]
    if not candidates:
        raise FileNotFoundError(f"No complete cache found under: {root}")

    def created_key(path: Path) -> tuple[str, str]:
        metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
        return str(metadata.get("created_utc", "")), path.name

    candidates.sort(key=created_key)
    chosen = candidates[-1]
    if len(candidates) > 1:
        print(f"Found {len(candidates)} complete caches; using newest: {chosen}")
    return chosen


def load_cache(cache_dir: Path):
    metadata = json.loads((cache_dir / "metadata.json").read_text(encoding="utf-8"))
    signature = metadata.get("signature", {})
    try:
        cpert_labels = [str(v) for v in signature["cpert"]]
        cact_labels = [str(v) for v in signature["cact"]]
    except KeyError as exc:
        raise ValueError(f"Cache metadata lacks parameter labels: {exc}") from exc

    s_raw = np.load(cache_dir / "s_raw.npy", mmap_mode="r")
    j_baseline = np.load(cache_dir / "j_baseline_mm.npy", mmap_mode="r")
    j_controlled = np.load(cache_dir / "j_controlled_mm.npy", mmap_mode="r")
    if s_raw.ndim != 4:
        raise ValueError(f"Expected s_raw[P,N,Y,X], got {s_raw.shape}")
    p_count, members, _, _ = s_raw.shape
    a_count = len(cact_labels)
    if len(cpert_labels) != p_count:
        raise ValueError("Cpert labels do not match s_raw")
    if j_baseline.shape != (p_count, members):
        raise ValueError(
            f"Expected j_baseline shape {(p_count, members)}, got {j_baseline.shape}"
        )
    if j_controlled.shape != (p_count, a_count, members):
        raise ValueError(
            "Expected j_controlled shape "
            f"{(p_count, a_count, members)}, got {j_controlled.shape}"
        )
    if not np.all(np.isfinite(j_baseline)) or not np.all(np.isfinite(j_controlled)):
        raise ValueError("Objective arrays contain non-finite values")
    return metadata, cpert_labels, cact_labels, s_raw, j_baseline, j_controlled


def nearest_index(values: np.ndarray, requested: float, name: str) -> int:
    index = int(np.argmin(np.abs(values - requested)))
    if not math.isclose(float(values[index]), float(requested), rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(f"Requested {name}={requested:g} is not present")
    return index


def subsample_weights(rng: np.random.Generator, replicates: int, members: int) -> np.ndarray:
    """Each row averages exactly 50 distinct members; subsets may overlap between rows."""
    if members < SAMPLE_SIZE:
        raise ValueError(f"Need at least {SAMPLE_SIZE} members, got {members}")
    weights = np.zeros((replicates, members), dtype=np.float64)
    for row in weights:
        row[rng.choice(members, size=SAMPLE_SIZE, replace=False)] = 1.0 / SAMPLE_SIZE
    return weights


def subset_weights(
    rng: np.random.Generator, repeats: int, members: int, subset_size: int
) -> np.ndarray:
    if subset_size == members:
        return np.full((1, members), 1.0 / members, dtype=np.float32)
    # Independent uniform subsets without replacement, one random ordering per row.
    scores = rng.random((repeats, members), dtype=np.float32)
    chosen = np.argpartition(scores, subset_size - 1, axis=1)[:, :subset_size]
    weights = np.zeros((repeats, members), dtype=np.float32)
    rows = np.repeat(np.arange(repeats, dtype=np.int32), subset_size)
    weights[rows, chosen.reshape(-1)] = np.float32(1.0 / subset_size)
    return weights


def vector_metrics(vectors: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference = np.asarray(reference, dtype=np.float64)
    reference_norm = float(np.linalg.norm(reference))
    if not np.isfinite(reference_norm) or reference_norm == 0.0:
        raise ValueError("Encountered a zero or non-finite reference-gradient norm")
    vectors64 = np.asarray(vectors, dtype=np.float64)
    vector_norms = np.linalg.norm(vectors64, axis=1)
    cosines = (vectors64 @ reference) / (vector_norms * reference_norm)
    angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    relative_l2 = np.linalg.norm(vectors64 - reference, axis=1) / reference_norm
    norm_ratio = vector_norms / reference_norm
    return angles, relative_l2, norm_ratio


def weighted_gradient_metrics(
    weights: np.ndarray,
    member_gradients: np.ndarray,
    reference: np.ndarray,
    deterministic: np.ndarray,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    count = weights.shape[0]
    angle = np.empty(count, dtype=np.float32)
    relative = np.empty(count, dtype=np.float32)
    norm_ratio = np.empty(count, dtype=np.float32)
    deflection = np.empty(count, dtype=np.float32)
    for start in range(0, count, batch_size):
        stop = min(start + batch_size, count)
        means = weights[start:stop] @ member_gradients
        a, r, n = vector_metrics(means, reference)
        d, _, _ = vector_metrics(means, deterministic)
        angle[start:stop] = a
        relative[start:stop] = r
        norm_ratio[start:stop] = n
        deflection[start:stop] = d
    return angle, relative, norm_ratio, deflection


def quantiles(values: np.ndarray) -> dict[str, float]:
    q = np.quantile(np.asarray(values, dtype=np.float64), [0.025, 0.5, 0.95, 0.975])
    return {
        "q025": float(q[0]),
        "median": float(q[1]),
        "q950": float(q[2]),
        "q975": float(q[3]),
    }


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_analysis(args: argparse.Namespace, cache_dir: Path, loaded) -> dict:
    metadata, cpert_labels, cact_labels, s_raw, j_baseline, j_controlled = loaded
    cpert = np.asarray([float(v) for v in cpert_labels], dtype=np.float64)
    cact = np.asarray([float(v) for v in cact_labels], dtype=np.float64)
    p_count, members, ny, nx = s_raw.shape
    if max(args.subset_sizes) > members:
        raise ValueError(f"Subset size exceeds available member count N={members}")
    subset_sizes = sorted(set(args.subset_sizes + [members]))
    deterministic_index = nearest_index(cpert, 0.0, "Cpert")
    highlight_index = nearest_index(cpert, args.highlight_cpert, "highlight Cpert")
    selected_cact_indices = [nearest_index(cact, v, "selected Cact") for v in args.selected_cact]

    seed_sequence = np.random.SeedSequence(args.seed)
    rng_subsample, rng_subset = [np.random.default_rng(s) for s in seed_sequence.spawn(2)]
    weights = subsample_weights(rng_subsample, args.subsample_replicates, members)

    # Positive accumulated precipitation reduction is beneficial: J_baseline - J_controlled.
    reduction = np.asarray(j_baseline[:, None, :] - j_controlled, dtype=np.float64)
    reduction_mean = reduction.mean(axis=2)
    reduction_member_sd = reduction.std(axis=2, ddof=1)
    response_subsample = np.einsum("bn,pan->pab", weights, reduction, optimize=True)
    response_ci_low = np.quantile(response_subsample, 0.025, axis=2)
    response_ci_high = np.quantile(response_subsample, 0.975, axis=2)
    response_probability_positive = (response_subsample > 0.0).mean(axis=2)
    best_index_subsample = np.argmax(response_subsample, axis=0)  # [Cact, subsample]
    best_probability = np.empty((p_count, len(cact)), dtype=np.float64)
    for p_index in range(p_count):
        best_probability[p_index] = (best_index_subsample == p_index).mean(axis=1)

    difference_vs_det = response_subsample - response_subsample[deterministic_index][None, :, :]
    difference_vs_det_ci_low = np.quantile(difference_vs_det, 0.025, axis=2)
    difference_vs_det_ci_high = np.quantile(difference_vs_det, 0.975, axis=2)
    probability_better_det = (difference_vs_det > 0.0).mean(axis=2)

    winner_rows = []
    for a_index, a_label in enumerate(cact_labels):
        winner = int(np.argmax(reduction_mean[:, a_index]))
        alternatives = np.delete(response_subsample[:, a_index, :], winner, axis=0)
        margin = response_subsample[winner, a_index] - alternatives.max(axis=0)
        q = quantiles(margin)
        winner_rows.append(
            {
                "cact": a_label,
                "full_sample_winner_cpert": cpert_labels[winner],
                "full_sample_winner_reduction_mm": reduction_mean[winner, a_index],
                "probability_full_sample_winner_remains_best": float(
                    (best_index_subsample[a_index] == winner).mean()
                ),
                "winner_margin_median_mm": q["median"],
                "winner_margin_ci_low_mm": q["q025"],
                "winner_margin_ci_high_mm": q["q975"],
            }
        )

    response_rows = []
    best_rows = []
    for p_index, p_label in enumerate(cpert_labels):
        for a_index, a_label in enumerate(cact_labels):
            response_rows.append(
                {
                    "cpert": p_label,
                    "cact": a_label,
                    "members": members,
                    "sample_size": SAMPLE_SIZE,
                    "full_sample_reduction_mean_mm": reduction_mean[p_index, a_index],
                    "subsample_reduction_median_mm": float(np.median(response_subsample[p_index, a_index])),
                    "member_sd_mm": reduction_member_sd[p_index, a_index],
                    "full_sample_standard_error_mm": reduction_member_sd[p_index, a_index] / math.sqrt(members),
                    "subsample_interval_low_mm": response_ci_low[p_index, a_index],
                    "subsample_interval_high_mm": response_ci_high[p_index, a_index],
                    "probability_positive_reduction": response_probability_positive[p_index, a_index],
                    "probability_best_cpert": best_probability[p_index, a_index],
                    "difference_vs_deterministic_mm": (
                        reduction_mean[p_index, a_index]
                        - reduction_mean[deterministic_index, a_index]
                    ),
                    "difference_vs_deterministic_ci_low_mm": difference_vs_det_ci_low[
                        p_index, a_index
                    ],
                    "difference_vs_deterministic_ci_high_mm": difference_vs_det_ci_high[
                        p_index, a_index
                    ],
                    "probability_better_than_deterministic": probability_better_det[
                        p_index, a_index
                    ],
                }
            )
            best_rows.append(
                {
                    "cact": a_label,
                    "cpert": p_label,
                    "probability_best": best_probability[p_index, a_index],
                }
            )

    print("Computing 50-member subset gradient variability...")
    gradients = np.asarray(s_raw.reshape(p_count, members, ny * nx), dtype=np.float32)
    full_gradient = gradients.mean(axis=1, dtype=np.float64)
    deterministic_gradient = full_gradient[deterministic_index]
    gradient_angle = np.empty((p_count, args.subsample_replicates), dtype=np.float32)
    gradient_relative = np.empty_like(gradient_angle)
    gradient_norm_ratio = np.empty_like(gradient_angle)
    gradient_deflection = np.empty_like(gradient_angle)
    gradient_subsample_rows = []
    for p_index, p_label in enumerate(cpert_labels):
        angle, relative, norm_ratio, deflection = weighted_gradient_metrics(
            weights,
            gradients[p_index],
            full_gradient[p_index],
            deterministic_gradient,
            args.batch_size,
        )
        gradient_angle[p_index] = angle
        gradient_relative[p_index] = relative
        gradient_norm_ratio[p_index] = norm_ratio
        gradient_deflection[p_index] = deflection
        qa, qr, qn, qd = map(quantiles, (angle, relative, norm_ratio, deflection))
        gradient_subsample_rows.append(
            {
                "cpert": p_label,
                "members": members,
                    "sample_size": SAMPLE_SIZE,
                "subsample_replicates": args.subsample_replicates,
                "angle_median_deg": qa["median"],
                "angle_ci_low_deg": qa["q025"],
                "angle_ci_high_deg": qa["q975"],
                "angle_95pct_deg": qa["q950"],
                "relative_l2_median": qr["median"],
                "relative_l2_ci_low": qr["q025"],
                "relative_l2_ci_high": qr["q975"],
                "relative_l2_95pct": qr["q950"],
                "norm_ratio_median": qn["median"],
                "norm_ratio_ci_low": qn["q025"],
                "norm_ratio_ci_high": qn["q975"],
                "deflection_from_deterministic_median_deg": qd["median"],
                "deflection_from_deterministic_ci_low_deg": qd["q025"],
                "deflection_from_deterministic_ci_high_deg": qd["q975"],
            }
        )

    print("Computing random-subset convergence...")
    subset_angle = np.full(
        (p_count, len(subset_sizes), args.subset_repeats), np.nan, dtype=np.float32
    )
    subset_relative = np.full_like(subset_angle, np.nan)
    subset_rows = []
    for size_index, size in enumerate(subset_sizes):
        size_weights = subset_weights(rng_subset, args.subset_repeats, members, size)
        actual_repeats = size_weights.shape[0]
        for p_index, p_label in enumerate(cpert_labels):
            angle, relative, _, _ = weighted_gradient_metrics(
                size_weights,
                gradients[p_index],
                full_gradient[p_index],
                deterministic_gradient,
                args.batch_size,
            )
            subset_angle[p_index, size_index, :actual_repeats] = angle
            subset_relative[p_index, size_index, :actual_repeats] = relative
            qa, qr = quantiles(angle), quantiles(relative)
            subset_rows.append(
                {
                    "cpert": p_label,
                    "subset_size": size,
                    "repeats": actual_repeats,
                    "angle_median_deg": qa["median"],
                    "angle_ci_low_deg": qa["q025"],
                    "angle_ci_high_deg": qa["q975"],
                    "angle_95pct_deg": qa["q950"],
                    "relative_l2_median": qr["median"],
                    "relative_l2_ci_low": qr["q025"],
                    "relative_l2_ci_high": qr["q975"],
                    "relative_l2_95pct": qr["q950"],
                }
            )

    print("Computing leave-one-out influence diagnostics...")
    gradient_jackknife_rows = []
    gradient_jackknife_max_angle = np.empty(p_count, dtype=np.float64)
    for p_index, p_label in enumerate(cpert_labels):
        leave_one_out = (
            members * full_gradient[p_index][None, :] - gradients[p_index]
        ) / (members - 1)
        angle, relative, norm_ratio = vector_metrics(leave_one_out, full_gradient[p_index])
        gradient_jackknife_max_angle[p_index] = float(angle.max())
        for member_index in range(members):
            gradient_jackknife_rows.append(
                {
                    "cpert": p_label,
                    "excluded_member": member_index + 1,
                    "angle_to_full_mean_deg": angle[member_index],
                    "relative_l2_difference": relative[member_index],
                    "norm_ratio": norm_ratio[member_index],
                }
            )

    response_jackknife_rows = []
    response_jackknife_max_shift = np.empty((p_count, len(cact)), dtype=np.float64)
    for p_index, p_label in enumerate(cpert_labels):
        for a_index, a_label in enumerate(cact_labels):
            leave_one_out = (
                members * reduction_mean[p_index, a_index] - reduction[p_index, a_index]
            ) / (members - 1)
            shifts = leave_one_out - reduction_mean[p_index, a_index]
            influence_member = int(np.argmax(np.abs(shifts)))
            response_jackknife_max_shift[p_index, a_index] = abs(shifts[influence_member])
            response_jackknife_rows.append(
                {
                    "cpert": p_label,
                    "cact": a_label,
                    "full_mean_reduction_mm": reduction_mean[p_index, a_index],
                    "max_absolute_leave_one_out_shift_mm": abs(shifts[influence_member]),
                    "most_influential_excluded_member": influence_member + 1,
                    "leave_one_out_min_reduction_mm": leave_one_out.min(),
                    "leave_one_out_max_reduction_mm": leave_one_out.max(),
                }
            )

    return {
        "source_metadata": metadata,
        "cache_dir": cache_dir,
        "cpert_labels": cpert_labels,
        "cact_labels": cact_labels,
        "cpert": cpert,
        "cact": cact,
        "members": members,
        "sample_size": SAMPLE_SIZE,
        "subset_sizes": np.asarray(subset_sizes, dtype=np.int32),
        "selected_cact_indices": selected_cact_indices,
        "highlight_index": highlight_index,
        "deterministic_index": deterministic_index,
        "reduction_mean": reduction_mean,
        "response_ci_low": response_ci_low,
        "response_ci_high": response_ci_high,
        "best_probability": best_probability,
        "response_subsample": response_subsample.astype(np.float32),
        "sample_member_indices": np.asarray([np.flatnonzero(row) + 1 for row in weights], dtype=np.int32),
        "gradient_angle": gradient_angle,
        "gradient_relative": gradient_relative,
        "gradient_norm_ratio": gradient_norm_ratio,
        "gradient_deflection": gradient_deflection,
        "subset_angle": subset_angle,
        "subset_relative": subset_relative,
        "gradient_jackknife_max_angle": gradient_jackknife_max_angle,
        "response_jackknife_max_shift": response_jackknife_max_shift,
        "response_rows": response_rows,
        "best_rows": best_rows,
        "winner_rows": winner_rows,
        "gradient_subsample_rows": gradient_subsample_rows,
        "subset_rows": subset_rows,
        "gradient_jackknife_rows": gradient_jackknife_rows,
        "response_jackknife_rows": response_jackknife_rows,
    }


def write_tables(output_dir: Path, result: dict) -> None:
    write_csv(
        output_dir / "response_subsample_summary.csv",
        list(result["response_rows"][0]),
        result["response_rows"],
    )
    write_csv(
        output_dir / "cpert_best_probability.csv",
        list(result["best_rows"][0]),
        result["best_rows"],
    )
    write_csv(
        output_dir / "winner_robustness.csv",
        list(result["winner_rows"][0]),
        result["winner_rows"],
    )
    write_csv(
        output_dir / "gradient_subsample_summary.csv",
        list(result["gradient_subsample_rows"][0]),
        result["gradient_subsample_rows"],
    )
    write_csv(
        output_dir / "gradient_subset_convergence.csv",
        list(result["subset_rows"][0]),
        result["subset_rows"],
    )
    write_csv(
        output_dir / "gradient_jackknife_memberwise.csv",
        list(result["gradient_jackknife_rows"][0]),
        result["gradient_jackknife_rows"],
    )
    write_csv(
        output_dir / "response_jackknife_summary.csv",
        list(result["response_jackknife_rows"][0]),
        result["response_jackknife_rows"],
    )

    np.savez_compressed(
        output_dir / "sampling_uncertainty.npz",
        sample_member_indices=result["sample_member_indices"],
        sample_size=SAMPLE_SIZE,
        cpert=result["cpert"],
        cact=result["cact"],
        subset_sizes=result["subset_sizes"],
        accumulated_precipitation_reduction_mean_mm=result["reduction_mean"],
        response_subsample_interval_low_mm=result["response_ci_low"],
        response_subsample_interval_high_mm=result["response_ci_high"],
        cpert_best_probability=result["best_probability"],
        response_subsample_reduction_mm=result["response_subsample"],
        gradient_subsample_angle_deg=result["gradient_angle"],
        gradient_subsample_relative_l2=result["gradient_relative"],
        gradient_subsample_norm_ratio=result["gradient_norm_ratio"],
        gradient_subsample_deflection_from_deterministic_deg=result["gradient_deflection"],
        subset_angle_deg=result["subset_angle"],
        subset_relative_l2=result["subset_relative"],
        gradient_jackknife_max_angle_deg=result["gradient_jackknife_max_angle"],
        response_jackknife_max_shift_mm=result["response_jackknife_max_shift"],
    )


def plot_results(output_dir: Path, result: dict) -> None:
    # Use a writable Matplotlib config directory when the default is read-only.
    # A task-specific temporary directory avoids a slow/repeated fallback and never
    # writes into the user's home directory.
    mpl_config = Path(tempfile.gettempdir()) / "ensemble_sampling_uncertainty_mplconfig"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Helvetica",
            "mathtext.it": "Helvetica:italic",
            "mathtext.bf": "Helvetica:bold",
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
        }
    )
    labels = result["cpert_labels"]
    x = np.arange(len(labels))
    highlight = result["highlight_index"]

    for legacy_name in (
        "control_response_subsample_interval.png",
        "gradient_subsample_uncertainty.png",
        "gradient_subset_convergence.png",
        "jackknife_member_influence.png",
    ):
        (output_dir / legacy_name).unlink(missing_ok=True)

    # Nonlinear response with paired-subsample percentile intervals.
    selected = result["selected_cact_indices"]
    selected_low = result["response_ci_low"][:, selected]
    selected_high = result["response_ci_high"][:, selected]
    y_min = float(np.min(selected_low))
    y_max = float(np.max(selected_high))
    y_padding = 0.04 * (y_max - y_min)
    for a_index in selected:
        fig, axis = plt.subplots(figsize=(5.2, 4.4))
        mean = np.median(result["response_subsample"][:, a_index, :], axis=1)
        low = result["response_ci_low"][:, a_index]
        high = result["response_ci_high"][:, a_index]
        # Points show median 50-member means; bars show their 2.5--97.5 percentiles.
        axis.vlines(x, low, high, color="#9cb4c4", lw=1.8, zorder=1)
        axis.hlines(low, x - 0.08, x + 0.08, color="#9cb4c4", lw=1.8, zorder=1)
        axis.hlines(high, x - 0.08, x + 0.08, color="#9cb4c4", lw=1.8, zorder=1)
        axis.plot(x, mean, "o-", color="#3f6f8f", lw=2.0, ms=6.5, zorder=2)
        axis.scatter(
            [highlight], [mean[highlight]], s=110, facecolors="white", edgecolors="#c62828", lw=2.5, zorder=5
        )
        axis.axhline(0.0, color="0.45", lw=1.0)
        axis.set_title(rf"$C_{{act}}={result['cact'][a_index]:g}$", fontsize=20, fontweight="bold")
        axis.set_xlabel(r"$C_{pert}$", fontsize=20)
        axis.set_ylabel("Accumulated precipitation\nreduction (mm)", fontsize=18)
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=55, ha="right", rotation_mode="anchor")
        axis.tick_params(axis="both", labelsize=16)
        axis.set_ylim(y_min - y_padding, y_max + y_padding)
        axis.grid(alpha=0.2)
        fig.subplots_adjust(left=0.23, right=0.97, bottom=0.27, top=0.88)
        fig.savefig(
            output_dir / f"control_response_subsample_interval_Cact_{result['cact'][a_index]:g}.png"
        )
        plt.close(fig)

    # Probability that each Cpert is the best among tested values.
    fig, axis = plt.subplots(figsize=(8.0, 6.0))
    image = axis.imshow(result["best_probability"], origin="lower", aspect="auto", vmin=0, vmax=1, cmap="viridis")
    for p_index in range(len(labels)):
        for a_index in range(len(result["cact_labels"])):
            value = result["best_probability"][p_index, a_index]
            color = "white" if value < 0.45 else "black"
            axis.text(a_index, p_index, f"{value:.2f}", ha="center", va="center", fontsize=14, color=color)
    axis.add_patch(
        plt.Rectangle((-0.5, highlight - 0.5), len(result["cact_labels"]), 1, fill=False, edgecolor="#e53935", lw=2)
    )
    axis.set_xticks(np.arange(len(result["cact_labels"])))
    axis.set_xticklabels(result["cact_labels"])
    axis.set_yticks(x)
    axis.set_yticklabels(labels)
    axis.tick_params(axis="both", labelsize=18)
    axis.set_xlabel(r"$C_{act}$", fontsize=20)
    axis.set_ylabel(r"$C_{pert}$", fontsize=20)
    axis.set_title(r"Subset probability of being the best $C_{pert}$", fontsize=20, fontweight="bold")
    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label("Probability", fontsize=18)
    colorbar.ax.tick_params(labelsize=18)
    fig.tight_layout()
    fig.savefig(output_dir / "cpert_best_probability.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    # Gradient subsample uncertainty versus Cpert.
    for values, ylabel, title, filename in (
        (
            result["gradient_angle"],
            r"Angle to full-$N$ mean (degree)",
            "50-member subset variability of gradient direction",
            "gradient_subsample_angle.png",
        ),
        (
            result["gradient_relative"],
            r"Relative $L_2$ difference",
            "50-member subset variability of ensemble-mean gradient",
            "gradient_subsample_relative_l2.png",
        ),
    ):
        fig, axis = plt.subplots(figsize=(6.2, 5.0))
        median = np.median(values, axis=1)
        low = np.quantile(values, 0.025, axis=1)
        high = np.quantile(values, 0.975, axis=1)
        axis.errorbar(
            x,
            median,
            yerr=np.vstack((median - low, high - median)),
            fmt="o-",
            color="#3f6f8f",
            ecolor="#9cb4c4",
            linewidth=2.0,
            markersize=7.0,
            elinewidth=1.8,
            capsize=4,
        )
        axis.scatter([highlight], [median[highlight]], s=110, facecolors="white", edgecolors="#c62828", lw=2.5, zorder=5)
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
        axis.tick_params(axis="both", labelsize=15)
        axis.set_xlabel(r"$C_{pert}$", fontsize=19)
        axis.set_ylabel(ylabel, fontsize=18)
        axis.set_title(title, fontsize=18, fontweight="bold")
        axis.grid(alpha=0.2)
        fig.subplots_adjust(left=0.18, right=0.97, bottom=0.22, top=0.88)
        fig.savefig(output_dir / filename, bbox_inches="tight")
        plt.close(fig)

    # Random-subset convergence. Show the median across random subsets.
    cmap = plt.get_cmap("viridis")
    for values, ylabel, filename in (
        (result["subset_angle"], r"Median angle to full-$N$ mean (degree)", "gradient_subset_convergence_angle.png"),
        (result["subset_relative"], r"Median relative $L_2$ difference", "gradient_subset_convergence_relative_l2.png"),
    ):
        fig, axis = plt.subplots(figsize=(5.6, 4.4))
        for p_index, p_label in enumerate(labels):
            color = "#c62828" if p_index == highlight else cmap(p_index / max(1, len(labels) - 1))
            width = 2.5 if p_index == highlight else 1.0
            zorder = 5 if p_index == highlight else 2
            median = np.nanmedian(values[p_index], axis=1)
            axis.plot(result["subset_sizes"], median, "o-", color=color, lw=width, ms=3.5, label=p_label, zorder=zorder)
        axis.set_xlabel("Subset size")
        axis.set_ylabel(ylabel)
        axis.set_xticks(result["subset_sizes"])
        axis.grid(alpha=0.2)
        axis.set_title("Random-subset convergence of ensemble-mean gradient")
        axis.legend(title=r"$C_{pert}$", ncol=2, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        fig.tight_layout()
        fig.savefig(output_dir / filename, bbox_inches="tight")
        plt.close(fig)

    # Leave-one-out influence.
    fig, axis = plt.subplots(figsize=(5.2, 4.2))
    axis.plot(x, result["gradient_jackknife_max_angle"], "o-", color="#3f6f8f")
    axis.scatter(
        [highlight], [result["gradient_jackknife_max_angle"][highlight]], s=70,
        facecolors="white", edgecolors="#c62828", lw=2, zorder=5
    )
    axis.set_xticks(x)
    axis.set_xticklabels(labels, rotation=45, ha="right")
    axis.set_xlabel(r"$C_{pert}$")
    axis.set_ylabel("Maximum leave-one-out angle (degree)")
    axis.set_title("Leave-one-out influence on ensemble-mean gradient")
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "gradient_jackknife_member_influence.png", bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.2, 4.8))
    image = axis.imshow(
        result["response_jackknife_max_shift"], origin="lower", aspect="auto", cmap="magma"
    )
    axis.set_xticks(np.arange(len(result["cact_labels"])))
    axis.set_xticklabels(result["cact_labels"])
    axis.set_yticks(x)
    axis.set_yticklabels(labels)
    axis.set_xlabel(r"$C_{act}$")
    axis.set_ylabel(r"$C_{pert}$")
    axis.set_title("Leave-one-out influence on accumulated precipitation response")
    fig.colorbar(image, ax=axis, label="Maximum mean-response shift (mm)")
    fig.tight_layout()
    fig.savefig(output_dir / "response_jackknife_member_influence.png", bbox_inches="tight")
    plt.close(fig)


def self_check() -> None:
    weights = subsample_weights(np.random.default_rng(7), 100, 100)
    assert weights.shape == (100, 100)
    assert np.all(np.count_nonzero(weights, axis=1) == SAMPLE_SIZE)
    assert np.allclose(weights.sum(axis=1), 1)
    assert np.all(weights[weights > 0] == 1.0 / SAMPLE_SIZE)
    assert np.array_equal(weights, subsample_weights(np.random.default_rng(7), 100, 100))
    values = np.arange(100, dtype=float)
    explicit = np.array([values[row > 0].mean() for row in weights])
    assert np.allclose(weights @ values, explicit)
    baseline = values + 100
    controlled = values * 0.5
    assert np.allclose(weights @ (baseline - controlled), weights @ baseline - weights @ controlled)
    reference = np.array([1.0, 2.0])
    angles, relative, _ = vector_metrics(np.tile(reference, (3, 1)), reference)
    assert np.allclose(angles, 0, atol=1e-5) and np.allclose(relative, 0)
    print("Self-check passed: exactly 50 unique members, paired means and reproducibility.")


def analysis_signature(args: argparse.Namespace, cache_dir: Path) -> dict:
    return {
        "version": 1,
        "sampling": "without replacement",
        "sample_size": SAMPLE_SIZE,
        "source": str(cache_dir),
        "source_files": {
            name: [int((cache_dir / name).stat().st_size), int((cache_dir / name).stat().st_mtime_ns)]
            for name in REQUIRED_CACHE_FILES
        },
        "seed": args.seed,
        "replicates": args.subsample_replicates,
        "subset_repeats": args.subset_repeats,
        "subset_sizes": args.subset_sizes,
        "highlight_cpert": args.highlight_cpert,
        "selected_cact": args.selected_cact,
        "batch_size": args.batch_size,
    }


def save_analysis_cache(path: Path, result: dict, signature: dict) -> None:
    path.mkdir(parents=True, exist_ok=True)
    arrays = {key: value for key, value in result.items() if isinstance(value, np.ndarray)}
    details = {key: str(value) if isinstance(value, Path) else value
               for key, value in result.items() if key not in arrays}
    np.savez_compressed(path / "arrays.npz", **arrays)
    (path / "details.json").write_text(
        json.dumps({"signature": signature, "result": details}, indent=2) + "\n",
        encoding="utf-8",
    )


def load_analysis_cache(path: Path, signature: dict) -> dict | None:
    if not (path / "arrays.npz").is_file() or not (path / "details.json").is_file():
        return None
    details = json.loads((path / "details.json").read_text(encoding="utf-8"))
    if details["signature"] != signature:
        return None
    result = details["result"]
    with np.load(path / "arrays.npz", allow_pickle=False) as arrays:
        result.update({key: arrays[key] for key in arrays.files})
    result["cache_dir"] = Path(result["cache_dir"])
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.self_check:
        self_check()
        return 0
    cache_dir = resolve_cache(args)
    loaded = load_cache(cache_dir)
    metadata, cpert_labels, cact_labels, s_raw, j_baseline, j_controlled = loaded
    if s_raw.shape[1] < SAMPLE_SIZE:
        raise ValueError(f"Need at least {SAMPLE_SIZE} source members")
    output_dir = args.output_dir.expanduser().resolve()
    analysis_cache_root = args.analysis_cache_root.expanduser().resolve()
    # Restrict writable paths to this script's subsampling result/cache directories.
    if not output_dir.is_relative_to(DEFAULT_OUTPUT_DIR.resolve()):
        raise ValueError(f"--output-dir must be within {DEFAULT_OUTPUT_DIR}")
    if not analysis_cache_root.is_relative_to(DEFAULT_ANALYSIS_CACHE_ROOT.resolve()):
        raise ValueError(f"--analysis-cache-root must be within {DEFAULT_ANALYSIS_CACHE_ROOT}")
    for destination in (output_dir, analysis_cache_root):
        if destination.is_relative_to(cache_dir) or cache_dir.is_relative_to(destination):
            raise ValueError("Writable directories must be separate from the source cache")
    print(f"Input cache (read-only): {cache_dir}")
    print(f"Raw gradients: {s_raw.shape}; baseline J: {j_baseline.shape}; controlled J: {j_controlled.shape}")
    print(f"Plan: {args.subsample_replicates} paired subsets of 50 distinct members; "
          f"{args.subset_repeats} convergence subsets per size")
    signature = analysis_signature(args, cache_dir)
    key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:16]
    analysis_cache_dir = analysis_cache_root / f"cache_{key}"
    print(f"Subsampling analysis cache: {analysis_cache_dir}")
    print(f"Subsampling results: {output_dir}")
    if args.dry_run:
        print("Dry run complete; nothing was written.")
        return 0
    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError("Result directory already contains files; use --overwrite to refresh it.")

    result = load_analysis_cache(analysis_cache_dir, signature)
    if result is None:
        if args.from_cache:
            raise FileNotFoundError(f"No matching 50-member analysis cache: {analysis_cache_dir}")
        result = run_analysis(args, cache_dir, loaded)
        save_analysis_cache(analysis_cache_dir, result, signature)
    else:
        print(f"Using matching analysis cache: {analysis_cache_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tables(output_dir, result)
    if not args.no_plots:
        print("Creating PNG figures...")
        plot_results(output_dir, result)
    output_metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "50-member random subsets without replacement",
        "input_cache_read_only": str(cache_dir),
        "analysis_cache": str(analysis_cache_dir),
        "signature": signature,
        "members_available": result["members"],
        "sample_size": SAMPLE_SIZE,
        "interval_interpretation": "2.5th--97.5th percentiles of 50-member subset means; not a bootstrap confidence interval",
        "pairing": "Same 50 member indices for baseline, controlled, gradients, and all Cpert/Cact combinations",
        "interpretation_limit": "Conditional on existing members and fixed control fields; not independent new-seed WRF experiments",
        "wrf_run": False,
        "original_workflow_modified": False,
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(output_metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Completed. Outputs: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
