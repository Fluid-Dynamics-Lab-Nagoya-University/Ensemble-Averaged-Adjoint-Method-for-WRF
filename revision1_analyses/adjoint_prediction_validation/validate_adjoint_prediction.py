#!/usr/bin/env python3
"""Quantitatively validate adjoint predictions against nonlinear WRF responses.

Simulation files are opened read-only. Analysis caches and results are written to
separate directories; replacing existing results requires --overwrite.

Data convention
---------------
* AD sensitivities and baseline states: dataset-root ``wrfout/AD/``
* Controlled nonlinear states: dataset-root ``wrfout/NL/``
* Objective: normalized 3x3 Gaussian-weighted terminal RAINNC
* Actual actuation: controlled QVAPOR(t=0) - baseline QVAPOR(t=0)

The supplied model's A_QVAPOR convention is converted to the objective-gradient
convention with ``raw_to_objective_scale`` (default: -1/60). Input locations can
be overridden with WRF_AD_ROOT/WRF_NL_ROOT or --ad-root/--nl-root.

by Shan Jiang, FDL, Nagoya University
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

# Loaded lazily after --dry-run so inventory checks require no scientific packages.
np = None
Dataset = None


DATASET_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AD_ROOT = Path(os.environ.get("WRF_AD_ROOT", str(DATASET_ROOT / "wrfout" / "AD")))
DEFAULT_NL_ROOT = Path(os.environ.get("WRF_NL_ROOT", str(DATASET_ROOT / "wrfout" / "NL")))
DEFAULT_REFERENCE_DELTAP = DATASET_ROOT / "dats" / "deltaP_mat_absG_i0M_NOV0.02_ngEQsg_sg_G3x3.npy"
DEFAULT_OUTPUT_ROOT = Path(__file__).resolve().parent / "results"
DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent / "cache"
DEFAULT_RUN_NAME = "main"
CACHE_SCHEMA_VERSION = 1

DEFAULT_CPERT = ("0", "0.001", "0.01", "0.05", "0.1", "0.15", "0.2", "0.3", "0.4", "0.5")
DEFAULT_CACT = ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9")

KERNEL = None


def load_scientific_dependencies() -> None:
    """Load numerical/NetCDF packages only when array analysis is requested."""
    global np, Dataset, KERNEL
    if np is not None and Dataset is not None and KERNEL is not None:
        return
    try:
        import numpy as numpy_module
    except ModuleNotFoundError as exc:  # pragma: no cover - runtime environment
        raise SystemExit(
            "Missing dependency 'numpy'. Run this script in the project's scientific "
            "Python environment; no package installation is performed automatically."
        ) from exc
    try:
        from netCDF4 import Dataset as dataset_class
    except ModuleNotFoundError as exc:  # pragma: no cover - runtime environment
        raise SystemExit(
            "Missing dependency 'netCDF4'. Run this script in the project's scientific "
            "Python environment; no package installation is performed automatically."
        ) from exc

    np = numpy_module
    Dataset = dataset_class
    KERNEL = np.array(
        [
            [0.25, 0.50, 0.25],
            [0.50, 1.00, 0.50],
            [0.25, 0.50, 0.25],
        ],
        dtype=np.float64,
    ) / 4.0


def canonical_decimal(value: str) -> str:
    """Return the directory-name form of a nonnegative decimal parameter."""
    try:
        number = Decimal(str(value).strip())
    except InvalidOperation as exc:
        raise argparse.ArgumentTypeError(f"Invalid decimal value: {value!r}") from exc
    if not number.is_finite() or number < 0:
        raise argparse.ArgumentTypeError(f"Expected a finite nonnegative value, got {value!r}")
    if number == 0:
        return "0"
    return format(number.normalize(), "f")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare ensemble-averaged first-order adjoint predictions with paired "
            "nonlinear WRF objective-function changes."
        )
    )
    parser.add_argument("--ad-root", type=Path, default=DEFAULT_AD_ROOT)
    parser.add_argument("--nl-root", type=Path, default=DEFAULT_NL_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument(
        "--reference-deltap",
        type=Path,
        default=DEFAULT_REFERENCE_DELTAP,
        help="Supplied dats/ 10x9 accumulated-precipitation-reduction matrix used only as a consistency check.",
    )
    parser.add_argument("--cpert", nargs="+", type=canonical_decimal, default=list(DEFAULT_CPERT))
    parser.add_argument("--cact", nargs="+", type=canonical_decimal, default=list(DEFAULT_CACT))
    parser.add_argument("--members", type=int, default=100)
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Parallel NetCDF reader processes used while building a cache (default: 4).",
    )
    parser.add_argument("--target-row", type=int, default=34, help="Zero-based first spatial index.")
    parser.add_argument("--target-col", type=int, default=16, help="Zero-based second spatial index.")
    parser.add_argument("--rain-time-index", type=int, default=6)
    parser.add_argument("--adjoint-time-index", type=int, default=-1)
    parser.add_argument(
        "--raw-to-objective-scale",
        type=float,
        default=-1.0 / 60.0,
        help="Convert raw A_QVAPOR inner products to the normalized 3x3 precipitation objective convention.",
    )
    parser.add_argument(
        "--error-floor-mm",
        type=float,
        default=0.1,
        help="Positive denominator floor used by robust normalized-error metrics.",
    )
    parser.add_argument(
        "--zero-tolerance",
        type=float,
        default=1.0e-12,
        help="Tolerance for the post-control QVAPOR zero-limit diagnostic.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help=(
            "Optional results-directory name (default: main). Existing results are "
            "never replaced unless --overwrite is also given."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow tables, metadata, NPZ, and figures in the selected result directory to be replaced.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Check all expected input paths without reading arrays or writing results.",
    )
    parser.add_argument(
        "--from-cache",
        action="store_true",
        help="Require and use an existing matching cache; never open the WRF NetCDF inputs.",
    )
    parser.add_argument("--no-plots", action="store_true", help="Write tables and NPZ only.")
    args = parser.parse_args(argv)

    if args.members <= 0:
        parser.error("--members must be positive")
    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.error_floor_mm <= 0:
        parser.error("--error-floor-mm must be positive")
    if args.zero_tolerance < 0:
        parser.error("--zero-tolerance must be nonnegative")
    if not math.isfinite(args.raw_to_objective_scale) or args.raw_to_objective_scale == 0:
        parser.error("--raw-to-objective-scale must be finite and nonzero")
    if args.run_name is not None:
        candidate = Path(args.run_name)
        if candidate.name != args.run_name or args.run_name in (".", ".."):
            parser.error("--run-name must be a single new directory name, not a path")
    if len(set(args.cpert)) != len(args.cpert):
        parser.error("--cpert contains duplicate values")
    if len(set(args.cact)) != len(args.cact):
        parser.error("--cact contains duplicate values")
    return args


def ad_directory(ad_root: Path, cpert: str) -> Path:
    return ad_root / f"woinput_absG_i0M_NOV0.02_sg_{cpert}"


def ad_file(ad_root: Path, cpert: str, member: int) -> Path:
    return ad_directory(ad_root, cpert) / (
        f"wrfout_d01_2018-07-05_120000_woinput_AD{member}"
    )


def nl_directory(nl_root: Path, cpert: str, cact: str) -> Path:
    return nl_root / f"absG_i0M_NOV0.02_ng{cpert}_ig{cact}_sg{cpert}"


def nl_file(nl_root: Path, cpert: str, cact: str, member: int) -> Path:
    return nl_directory(nl_root, cpert, cact) / (
        f"wrfout_d01_2018-07-05_120000_woinput_NL{member}"
    )


def validate_input_paths(args: argparse.Namespace) -> None:
    """Strictly check paired AD/NL names with one directory listing per case."""
    if not args.ad_root.is_dir():
        raise FileNotFoundError(f"AD root does not exist: {args.ad_root}")
    if not args.nl_root.is_dir():
        raise FileNotFoundError(f"NL root does not exist: {args.nl_root}")

    missing: List[Path] = []

    def check_directory(directory: Path, expected_names: Sequence[str]) -> None:
        try:
            available_names = set(os.listdir(directory))
        except FileNotFoundError:
            missing.extend(directory / name for name in expected_names)
            return
        missing.extend(
            directory / name for name in expected_names if name not in available_names
        )

    for cpert in args.cpert:
        print(f"  inventory: Cpert={cpert}", flush=True)
        check_directory(
            ad_directory(args.ad_root, cpert),
            [
                f"wrfout_d01_2018-07-05_120000_woinput_AD{member}"
                for member in range(1, args.members + 1)
            ],
        )
        for cact in args.cact:
            check_directory(
                nl_directory(args.nl_root, cpert, cact),
                [
                    f"wrfout_d01_2018-07-05_120000_woinput_NL{member}"
                    for member in range(1, args.members + 1)
                ],
            )

    if missing:
        preview = "\n".join(f"  {path}" for path in missing[:20])
        remainder = len(missing) - min(20, len(missing))
        suffix = f"\n  ... and {remainder} more" if remainder else ""
        raise FileNotFoundError(
            f"Missing {len(missing)} required paired input files:\n{preview}{suffix}"
        )


def require_variable(dataset: Dataset, name: str, path: Path):
    if name not in dataset.variables:
        raise KeyError(f"Variable {name!r} is missing from {path}")
    return dataset.variables[name]


def read_surface_field(variable, time_index: int, variable_name: str, path: Path) -> np.ndarray:
    if variable.ndim != 4:
        raise ValueError(
            f"Expected 4-D {variable_name} in {path}, got shape {variable.shape}"
        )
    try:
        field = variable[time_index, 0, :, :]
    except IndexError as exc:
        raise IndexError(
            f"Time index {time_index} is invalid for {variable_name} shape {variable.shape} in {path}"
        ) from exc
    return np.asarray(field, dtype=np.float64)


def read_rain_field(variable, time_index: int, path: Path) -> np.ndarray:
    if variable.ndim != 3:
        raise ValueError(f"Expected 3-D RAINNC in {path}, got shape {variable.shape}")
    try:
        field = variable[time_index, :, :]
    except IndexError as exc:
        raise IndexError(
            f"Time index {time_index} is invalid for RAINNC shape {variable.shape} in {path}"
        ) from exc
    return np.asarray(field, dtype=np.float64)


def objective(rain: np.ndarray, target_row: int, target_col: int) -> float:
    """Evaluate the normalized 3x3 weighted RAINNC objective at the target."""
    if rain.ndim != 2:
        raise ValueError(f"Expected a 2-D RAINNC field, got shape {rain.shape}")
    row0, row1 = target_row - 1, target_row + 2
    col0, col1 = target_col - 1, target_col + 2
    if row0 < 0 or col0 < 0 or row1 > rain.shape[0] or col1 > rain.shape[1]:
        raise ValueError(
            f"Target ({target_row}, {target_col}) is invalid for RAINNC shape {rain.shape}"
        )
    return float(np.sum(rain[row0:row1, col0:col1] * KERNEL))


def read_ad_member(
    path: Path,
    rain_time_index: int,
    adjoint_time_index: int,
    target_row: int,
    target_col: int,
) -> Dict[str, object]:
    with Dataset(path, mode="r") as dataset:
        sensitivity = read_surface_field(
            require_variable(dataset, "A_QVAPOR", path),
            adjoint_time_index,
            "A_QVAPOR",
            path,
        )
        qvapor = read_surface_field(
            require_variable(dataset, "QVAPOR", path), 0, "QVAPOR", path
        )
        rain = read_rain_field(
            require_variable(dataset, "RAINNC", path), rain_time_index, path
        )
    if sensitivity.shape != qvapor.shape:
        raise ValueError(
            f"A_QVAPOR/QVAPOR shape mismatch in {path}: {sensitivity.shape} vs {qvapor.shape}"
        )
    return {
        "s_raw": sensitivity,
        "q_baseline": qvapor,
        "j_baseline": objective(rain, target_row, target_col),
    }


def read_nl_member(
    path: Path,
    rain_time_index: int,
    target_row: int,
    target_col: int,
) -> Dict[str, object]:
    with Dataset(path, mode="r") as dataset:
        qvapor = read_surface_field(
            require_variable(dataset, "QVAPOR", path), 0, "QVAPOR", path
        )
        rain = read_rain_field(
            require_variable(dataset, "RAINNC", path), rain_time_index, path
        )
    return {
        "q_controlled": qvapor,
        "j_controlled": objective(rain, target_row, target_col),
    }


def error_metrics(predicted: float, true: float, floor: float) -> Dict[str, float]:
    absolute = abs(true - predicted)
    return {
        "absolute_error_mm": absolute,
        "truth_normalized_error": absolute / max(abs(true), floor),
        "symmetric_normalized_error": (
            2.0 * absolute / max(abs(true) + abs(predicted), floor)
        ),
    }


def sample_std(values: Sequence[float]) -> float:
    return float(np.std(np.asarray(values, dtype=np.float64), ddof=1)) if len(values) > 1 else 0.0


def read_ad_task(task: Tuple[object, ...]) -> Tuple[int, np.ndarray, np.ndarray, float]:
    """Process-pool worker: read one AD member and return only required slices."""
    load_scientific_dependencies()
    member, path_string, rain_index, adjoint_index, target_row, target_col = task
    data = read_ad_member(
        Path(path_string), rain_index, adjoint_index, target_row, target_col
    )
    return (
        int(member),
        np.asarray(data["s_raw"], dtype=np.float32),
        np.asarray(data["q_baseline"], dtype=np.float32),
        float(data["j_baseline"]),
    )


def read_nl_task(task: Tuple[object, ...]) -> Tuple[int, np.ndarray, float]:
    """Process-pool worker: read one nonlinear member and return required slices."""
    load_scientific_dependencies()
    member, path_string, rain_index, target_row, target_col = task
    data = read_nl_member(Path(path_string), rain_index, target_row, target_col)
    return (
        int(member),
        np.asarray(data["q_controlled"], dtype=np.float32),
        float(data["j_controlled"]),
    )


def run_reader_tasks(executor, worker, tasks: Sequence[Tuple[object, ...]]):
    if executor is None:
        return [worker(task) for task in tasks]
    return list(executor.map(worker, tasks, chunksize=1))


def cache_signature_payload(args: argparse.Namespace) -> Dict[str, object]:
    """Configuration fields that determine the extracted cache contents."""
    return {
        "schema_version": CACHE_SCHEMA_VERSION,
        "ad_root": str(args.ad_root.resolve()),
        "nl_root": str(args.nl_root.resolve()),
        "cpert": list(args.cpert),
        "cact": list(args.cact),
        "members": args.members,
        "target_row": args.target_row,
        "target_col": args.target_col,
        "rain_time_index": args.rain_time_index,
        "adjoint_time_index": args.adjoint_time_index,
        "zero_tolerance": args.zero_tolerance,
    }


def matching_cache_directory(args: argparse.Namespace) -> Path:
    payload = cache_signature_payload(args)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:16]
    return args.cache_root / f"cache_{digest}"


def extract_cache_arrays(args: argparse.Namespace) -> Dict[str, np.ndarray]:
    """Read required NetCDF slices in parallel and retain compact arrays in memory."""
    s_raw_by_cpert: List[np.ndarray] = []
    j_baseline_by_cpert: List[np.ndarray] = []
    delta_by_cpert: List[np.ndarray] = []
    j_controlled_by_cpert: List[np.ndarray] = []
    new_zero_by_cpert: List[np.ndarray] = []

    executor = None
    if args.workers > 1:
        executor = ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=multiprocessing.get_context("spawn"),
        )
    try:
        for cpert in args.cpert:
            print(
                f"[cache] parallel AD read: Cpert={cpert}, members={args.members}, "
                f"workers={args.workers}",
                flush=True,
            )
            ad_tasks = [
                (
                    member,
                    str(ad_file(args.ad_root, cpert, member)),
                    args.rain_time_index,
                    args.adjoint_time_index,
                    args.target_row,
                    args.target_col,
                )
                for member in range(1, args.members + 1)
            ]
            ad_results = run_reader_tasks(executor, read_ad_task, ad_tasks)
            if [result[0] for result in ad_results] != list(range(1, args.members + 1)):
                raise RuntimeError(f"AD member ordering failed for Cpert={cpert}")
            s_raw = np.stack([result[1] for result in ad_results], axis=0)
            q_baseline = np.stack([result[2] for result in ad_results], axis=0)
            j_baseline = np.asarray([result[3] for result in ad_results], dtype=np.float64)

            delta_by_cact: List[np.ndarray] = []
            j_controlled_by_cact: List[np.ndarray] = []
            new_zero_by_cact: List[np.ndarray] = []
            for cact in args.cact:
                print(
                    f"  [cache] parallel NL read: Cpert={cpert}, Cact={cact}",
                    flush=True,
                )
                nl_tasks = [
                    (
                        member,
                        str(nl_file(args.nl_root, cpert, cact, member)),
                        args.rain_time_index,
                        args.target_row,
                        args.target_col,
                    )
                    for member in range(1, args.members + 1)
                ]
                nl_results = run_reader_tasks(executor, read_nl_task, nl_tasks)
                if [result[0] for result in nl_results] != list(range(1, args.members + 1)):
                    raise RuntimeError(
                        f"NL member ordering failed for Cpert={cpert}, Cact={cact}"
                    )
                q_controlled = np.stack([result[1] for result in nl_results], axis=0)
                if q_controlled.shape != q_baseline.shape:
                    raise ValueError(
                        f"Baseline/controlled QVAPOR shape mismatch for Cpert={cpert}, "
                        f"Cact={cact}: {q_baseline.shape} vs {q_controlled.shape}"
                    )
                delta_actual = q_controlled - q_baseline
                newly_zero = (
                    (q_controlled <= args.zero_tolerance)
                    & (q_baseline > args.zero_tolerance)
                    & (delta_actual < -args.zero_tolerance)
                )
                new_zero_fraction = np.count_nonzero(newly_zero, axis=(1, 2)) / float(
                    newly_zero.shape[1] * newly_zero.shape[2]
                )
                delta_by_cact.append(np.asarray(delta_actual, dtype=np.float32))
                j_controlled_by_cact.append(
                    np.asarray([result[2] for result in nl_results], dtype=np.float64)
                )
                new_zero_by_cact.append(
                    np.asarray(new_zero_fraction, dtype=np.float64)
                )

            s_raw_by_cpert.append(s_raw)
            j_baseline_by_cpert.append(j_baseline)
            delta_by_cpert.append(np.stack(delta_by_cact, axis=0))
            j_controlled_by_cpert.append(np.stack(j_controlled_by_cact, axis=0))
            new_zero_by_cpert.append(np.stack(new_zero_by_cact, axis=0))
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    return {
        "s_raw": np.stack(s_raw_by_cpert, axis=0),
        "j_baseline_mm": np.stack(j_baseline_by_cpert, axis=0),
        "delta_actual": np.stack(delta_by_cpert, axis=0),
        "j_controlled_mm": np.stack(j_controlled_by_cpert, axis=0),
        "new_zero_grid_fraction": np.stack(new_zero_by_cpert, axis=0),
    }


def write_npy_exclusive(path: Path, array: np.ndarray) -> None:
    with path.open("xb") as handle:
        np.save(handle, array, allow_pickle=False)


def write_cache(
    cache_dir: Path, args: argparse.Namespace, arrays: Mapping[str, np.ndarray]
) -> None:
    """Create a complete new cache directory; never replace an existing cache."""
    ensure_output_is_separate(args.cache_root, (args.ad_root, args.nl_root))
    args.cache_root.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=False, exist_ok=False)
    for name, array in arrays.items():
        write_npy_exclusive(cache_dir / f"{name}.npy", array)
    metadata = {
        "complete": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "signature": cache_signature_payload(args),
        "workers_used_for_extraction": args.workers,
        "arrays": {
            name: {"shape": list(array.shape), "dtype": str(array.dtype)}
            for name, array in arrays.items()
        },
        "input_files_opened_read_only": True,
        "existing_simulation_data_modified": False,
    }
    with (cache_dir / "metadata.json").open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def load_cache(cache_dir: Path, args: argparse.Namespace) -> Dict[str, np.ndarray]:
    metadata_path = cache_dir / "metadata.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(
            f"Cache is missing or incomplete (no metadata.json): {cache_dir}"
        )
    with metadata_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("complete") is not True:
        raise ValueError(f"Cache is not marked complete: {cache_dir}")
    expected_signature = cache_signature_payload(args)
    if metadata.get("signature") != expected_signature:
        raise ValueError(f"Cache signature does not match requested analysis: {cache_dir}")

    arrays: Dict[str, np.ndarray] = {}
    for name in (
        "s_raw",
        "j_baseline_mm",
        "delta_actual",
        "j_controlled_mm",
        "new_zero_grid_fraction",
    ):
        path = cache_dir / f"{name}.npy"
        if not path.is_file():
            raise FileNotFoundError(f"Cache array is missing: {path}")
        arrays[name] = np.load(path, mmap_mode="r", allow_pickle=False)
    return arrays


def obtain_cache(args: argparse.Namespace) -> Tuple[Dict[str, np.ndarray], Path, str]:
    cache_dir = matching_cache_directory(args)
    if cache_dir.exists():
        print(f"Using matching read-only analysis cache: {cache_dir}", flush=True)
        return load_cache(cache_dir, args), cache_dir, "reused"
    if args.from_cache:
        raise FileNotFoundError(
            f"--from-cache was requested, but no matching cache exists: {cache_dir}"
        )

    print(f"No matching cache found; extracting with {args.workers} workers.", flush=True)
    arrays = extract_cache_arrays(args)
    try:
        write_cache(cache_dir, args, arrays)
    except PermissionError:
        print(
            f"Permission denied while creating cache under {args.cache_root}; no input "
            "simulation data were modified.",
            file=sys.stderr,
        )
        raise
    print(f"New non-overwriting cache written to: {cache_dir}", flush=True)
    return arrays, cache_dir, "created"


def analyze_from_cache(
    args: argparse.Namespace, cache: Mapping[str, np.ndarray]
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    """Vectorized member calculations from compact in-memory or memory-mapped arrays."""
    s_raw_all = cache["s_raw"]
    j_baseline_all = cache["j_baseline_mm"]
    delta_all = cache["delta_actual"]
    j_controlled_all = cache["j_controlled_mm"]
    new_zero_all = cache["new_zero_grid_fraction"]
    expected_leading = (len(args.cpert), len(args.cact), args.members)
    if delta_all.shape[:3] != expected_leading:
        raise ValueError(
            f"Cache delta_actual shape {delta_all.shape} is incompatible with {expected_leading}"
        )

    member_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    for i_cpert, cpert in enumerate(args.cpert):
        s_raw = s_raw_all[i_cpert]
        s_raw_mean = np.mean(s_raw, axis=0, dtype=np.float64)
        j_baseline = np.asarray(j_baseline_all[i_cpert], dtype=np.float64)
        for i_cact, cact in enumerate(args.cact):
            delta_actual = delta_all[i_cpert, i_cact]
            j_controlled = np.asarray(
                j_controlled_all[i_cpert, i_cact], dtype=np.float64
            )
            new_zero_fraction = np.asarray(
                new_zero_all[i_cpert, i_cact], dtype=np.float64
            )
            member_specific_raw_dots = np.einsum(
                "kij,kij->k",
                s_raw,
                delta_actual,
                dtype=np.float64,
                casting="unsafe",
                optimize=True,
            )
            mean_sensitivity_raw_dots = np.einsum(
                "ij,kij->k",
                s_raw_mean,
                delta_actual,
                dtype=np.float64,
                casting="unsafe",
                optimize=True,
            )
            member_specific_predictions = (
                args.raw_to_objective_scale * member_specific_raw_dots
            )
            predictions = args.raw_to_objective_scale * mean_sensitivity_raw_dots
            true_changes = j_controlled - j_baseline
            absolute_errors = np.abs(true_changes - predictions)
            truth_errors = absolute_errors / np.maximum(
                np.abs(true_changes), args.error_floor_mm
            )
            symmetric_errors = 2.0 * absolute_errors / np.maximum(
                np.abs(true_changes) + np.abs(predictions), args.error_floor_mm
            )
            flat_delta = delta_actual.reshape(args.members, -1)
            actuation_l1 = np.sum(np.abs(flat_delta), axis=1, dtype=np.float64)
            actuation_l2 = np.linalg.norm(flat_delta.astype(np.float64), axis=1)
            actuation_linf = np.max(np.abs(flat_delta), axis=1)

            for member_index in range(args.members):
                member_rows.append(
                    {
                        "cpert": float(cpert),
                        "cpert_label": cpert,
                        "cact": float(cact),
                        "cact_label": cact,
                        "member": member_index + 1,
                        "j_baseline_mm": float(j_baseline[member_index]),
                        "j_controlled_mm": float(j_controlled[member_index]),
                        "raw_s_dot_delta": float(
                            mean_sensitivity_raw_dots[member_index]
                        ),
                        "raw_member_specific_s_dot_delta": float(
                            member_specific_raw_dots[member_index]
                        ),
                        "delta_j_pred_mm": float(predictions[member_index]),
                        "delta_j_pred_member_specific_mm": float(
                            member_specific_predictions[member_index]
                        ),
                        "delta_j_true_mm": float(true_changes[member_index]),
                        "rainfall_reduction_pred_mm": float(-predictions[member_index]),
                        "rainfall_reduction_true_mm": float(-true_changes[member_index]),
                        "actuation_l1": float(actuation_l1[member_index]),
                        "actuation_l2": float(actuation_l2[member_index]),
                        "actuation_linf": float(actuation_linf[member_index]),
                        "new_zero_grid_fraction": float(new_zero_fraction[member_index]),
                        "absolute_error_mm": float(absolute_errors[member_index]),
                        "truth_normalized_error": float(truth_errors[member_index]),
                        "symmetric_normalized_error": float(
                            symmetric_errors[member_index]
                        ),
                    }
                )

            delta_mean = np.mean(delta_actual, axis=0, dtype=np.float64)
            ensemble_mean_sensitivity_prediction = float(np.mean(predictions))
            paired_member_specific_prediction = float(
                np.mean(member_specific_predictions)
            )
            true_change = float(np.mean(true_changes))
            common_mean_prediction = float(
                args.raw_to_objective_scale
                * np.einsum(
                    "ij,ij->",
                    s_raw_mean,
                    delta_mean,
                    dtype=np.float64,
                    optimize=True,
                )
            )
            mean_delta_norm = float(np.linalg.norm(delta_mean.ravel()))
            centered_delta = flat_delta.astype(np.float64) - delta_mean.ravel()
            member_variability = float(
                np.mean(np.linalg.norm(centered_delta, axis=1))
                / max(float(np.mean(actuation_l2)), np.finfo(np.float64).tiny)
            )
            summary_rows.append(
                {
                    "cpert": float(cpert),
                    "cpert_label": cpert,
                    "cact": float(cact),
                    "cact_label": cact,
                    "n_members": args.members,
                    "delta_j_pred_mm": ensemble_mean_sensitivity_prediction,
                    "delta_j_true_mm": true_change,
                    "rainfall_reduction_pred_mm": -ensemble_mean_sensitivity_prediction,
                    "rainfall_reduction_true_mm": -true_change,
                    "delta_j_pred_member_std_mm": sample_std(predictions),
                    "delta_j_true_member_std_mm": sample_std(true_changes),
                    "delta_j_pred_common_mean_actuation_mm": common_mean_prediction,
                    "delta_j_pred_paired_member_specific_mm": (
                        paired_member_specific_prediction
                    ),
                    "paired_minus_common_prediction_mm": (
                        paired_member_specific_prediction
                        - ensemble_mean_sensitivity_prediction
                    ),
                    "mean_actual_actuation_l2": float(np.mean(actuation_l2)),
                    "mean_actual_actuation_field_l2": mean_delta_norm,
                    "actual_actuation_member_variability": member_variability,
                    "mean_new_zero_grid_fraction": float(np.mean(new_zero_fraction)),
                    "members_with_new_zero_fraction": float(
                        np.mean(new_zero_fraction > 0.0)
                    ),
                    **error_metrics(
                        ensemble_mean_sensitivity_prediction,
                        true_change,
                        args.error_floor_mm,
                    ),
                }
            )

    deterministic_rows: List[Dict[str, object]] = []
    if "0" in args.cpert:
        for cact in args.cact:
            source = next(
                row
                for row in member_rows
                if row["cpert_label"] == "0"
                and row["cact_label"] == cact
                and row["member"] == 1
            )
            deterministic_rows.append(
                {
                    key: source[key]
                    for key in (
                        "cact",
                        "cact_label",
                        "member",
                        "delta_j_pred_mm",
                        "delta_j_true_mm",
                        "rainfall_reduction_pred_mm",
                        "rainfall_reduction_true_mm",
                        "absolute_error_mm",
                        "truth_normalized_error",
                        "symmetric_normalized_error",
                        "actuation_l2",
                        "new_zero_grid_fraction",
                    )
                }
            )
    return member_rows, summary_rows, deterministic_rows


def add_reference_deltap(
    summary_rows: List[Dict[str, object]], reference_path: Path
) -> Dict[str, object]:
    report: Dict[str, object] = {
        "path": str(reference_path),
        "available": reference_path.is_file(),
    }
    if not reference_path.is_file():
        for row in summary_rows:
            row["existing_delta_p_mm"] = math.nan
            row["new_minus_existing_reduction_mm"] = math.nan
        return report

    matrix = np.load(reference_path)
    expected_shape = (len(DEFAULT_CPERT), len(DEFAULT_CACT))
    if matrix.shape != expected_shape:
        raise ValueError(
            f"Reference deltaP matrix {reference_path} has shape {matrix.shape}; "
            f"expected {expected_shape}"
        )
    cpert_index = {label: index for index, label in enumerate(DEFAULT_CPERT)}
    cact_index = {label: index for index, label in enumerate(DEFAULT_CACT)}
    differences: List[float] = []
    for row in summary_rows:
        cp = str(row["cpert_label"])
        ca = str(row["cact_label"])
        if cp in cpert_index and ca in cact_index:
            existing = float(matrix[cpert_index[cp], cact_index[ca]])
            difference = float(row["rainfall_reduction_true_mm"] - existing)
            differences.append(difference)
        else:
            existing = math.nan
            difference = math.nan
        row["existing_delta_p_mm"] = existing
        row["new_minus_existing_reduction_mm"] = difference

    report.update(
        {
            "comparisons": len(differences),
            "max_abs_difference_mm": (
                max(abs(value) for value in differences) if differences else math.nan
            ),
            "mean_abs_difference_mm": (
                float(np.mean(np.abs(differences))) if differences else math.nan
            ),
        }
    )
    return report


def linear_fit_metrics(rows: Sequence[Mapping[str, object]]) -> Dict[str, float]:
    predicted = np.asarray([row["delta_j_pred_mm"] for row in rows], dtype=np.float64)
    true = np.asarray([row["delta_j_true_mm"] for row in rows], dtype=np.float64)
    if len(rows) < 2:
        return {
            "n": float(len(rows)),
            "slope": math.nan,
            "intercept_mm": math.nan,
            "r_squared": math.nan,
            "pearson_r": math.nan,
            "mae_mm": float(np.mean(np.abs(true - predicted))) if len(rows) else math.nan,
            "rmse_mm": float(np.sqrt(np.mean((true - predicted) ** 2))) if len(rows) else math.nan,
        }
    slope, intercept = np.polyfit(predicted, true, 1)
    fitted = slope * predicted + intercept
    residual_sum = float(np.sum((true - fitted) ** 2))
    total_sum = float(np.sum((true - np.mean(true)) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0 else math.nan
    pearson_r = float(np.corrcoef(predicted, true)[0, 1])
    return {
        "n": float(len(rows)),
        "slope": float(slope),
        "intercept_mm": float(intercept),
        "r_squared": r_squared,
        "pearson_r": pearson_r,
        "mae_mm": float(np.mean(np.abs(true - predicted))),
        "rmse_mm": float(np.sqrt(np.mean((true - predicted) ** 2))),
    }


def ensure_output_is_separate(output_root: Path, input_roots: Sequence[Path]) -> None:
    output_resolved = output_root.expanduser().resolve()
    for input_root in input_roots:
        input_resolved = input_root.expanduser().resolve()
        try:
            common = Path(os.path.commonpath([str(output_resolved), str(input_resolved)]))
        except ValueError:
            continue
        if common == input_resolved:
            raise ValueError(
                f"Output root {output_resolved} must not be inside input-data root {input_resolved}"
            )


def new_run_directory(
    output_root: Path, run_name: str | None, overwrite: bool = False
) -> Tuple[Path, bool]:
    """Return the output directory and whether existing files may be replaced."""
    ensure_output_is_separate(output_root, (DEFAULT_AD_ROOT, DEFAULT_NL_ROOT))
    output_root.mkdir(parents=True, exist_ok=True)
    name = run_name or DEFAULT_RUN_NAME
    run_dir = output_root / name
    replacing_existing = run_dir.exists()
    if replacing_existing:
        if not run_dir.is_dir():
            raise FileExistsError(f"Result path exists and is not a directory: {run_dir}")
        if not overwrite:
            raise FileExistsError(
                f"Result directory already exists: {run_dir}; pass --overwrite to replace its outputs"
            )
    else:
        run_dir.mkdir(parents=False, exist_ok=False)
    return run_dir, replacing_existing


def write_csv(
    path: Path, rows: Sequence[Mapping[str, object]], overwrite: bool = False
) -> None:
    if not rows:
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    mode = "w" if overwrite else "x"
    with path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summary_arrays(
    summary_rows: Sequence[Mapping[str, object]],
    cpert: Sequence[str],
    cact: Sequence[str],
) -> Dict[str, np.ndarray]:
    keys = (
        "delta_j_pred_mm",
        "delta_j_true_mm",
        "rainfall_reduction_pred_mm",
        "rainfall_reduction_true_mm",
        "absolute_error_mm",
        "truth_normalized_error",
        "symmetric_normalized_error",
        "delta_j_pred_member_std_mm",
        "delta_j_true_member_std_mm",
        "delta_j_pred_common_mean_actuation_mm",
        "delta_j_pred_paired_member_specific_mm",
        "paired_minus_common_prediction_mm",
        "mean_new_zero_grid_fraction",
        "members_with_new_zero_fraction",
    )
    row_index = {label: index for index, label in enumerate(cpert)}
    col_index = {label: index for index, label in enumerate(cact)}
    arrays = {
        key: np.full((len(cpert), len(cact)), np.nan, dtype=np.float64) for key in keys
    }
    for row in summary_rows:
        i = row_index[str(row["cpert_label"])]
        j = col_index[str(row["cact_label"])]
        for key in keys:
            arrays[key][i, j] = float(row[key])
    arrays["cpert"] = np.asarray([float(value) for value in cpert], dtype=np.float64)
    arrays["cact"] = np.asarray([float(value) for value in cact], dtype=np.float64)
    return arrays


def save_npz(
    path: Path,
    summary_rows: Sequence[Mapping[str, object]],
    cpert: Sequence[str],
    cact: Sequence[str],
    overwrite: bool = False,
) -> None:
    arrays = summary_arrays(summary_rows, cpert, cact)
    mode = "wb" if overwrite else "xb"
    with path.open(mode) as handle:
        np.savez_compressed(handle, **arrays)


def plotting_module():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:  # pragma: no cover - runtime environment
        raise SystemExit(
            "Missing dependency 'matplotlib'. Use --no-plots or run in the project's "
            "scientific Python environment."
        ) from exc
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica"],
            "font.size": 18,
            "axes.labelsize": 22,
            "axes.titlesize": 20,
            "axes.linewidth": 1.5,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "legend.fontsize": 16,
            "lines.linewidth": 2.2,
            "lines.markersize": 7,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Helvetica",
            "mathtext.it": "Helvetica:italic",
            "mathtext.bf": "Helvetica:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    return plt


def new_plot_canvas(plt):
    """Create the same canvas and square plotting frame for every figure."""
    fig = plt.figure(figsize=(6.2, 8.5))
    ax = fig.add_axes((0.17, 0.47, 0.80, 0.48))
    ax.set_box_aspect(1)
    return fig, ax


def add_figure_legend(fig, ax) -> None:
    """Place every legend in the same reserved area below the plotting frame."""
    handles, labels = ax.get_legend_handles_labels()
    ncol = 3 if len(handles) > 12 else 2
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.33),
        ncol=ncol,
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.5,
        borderaxespad=0.0,
    )


def save_figure(fig, run_dir: Path, stem: str, overwrite: bool = False) -> None:
    path = run_dir / f"{stem}.png"
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing figure: {path}")
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=0.05)


def rotate_overlapping_xticklabels(fig, ax, angle: float = 45) -> bool:
    """Rotate x tick labels only when adjacent rendered labels overlap."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = [label for label in ax.get_xticklabels() if label.get_visible() and label.get_text()]
    boxes = [label.get_window_extent(renderer) for label in labels]
    if not any(left.x1 > right.x0 for left, right in zip(boxes, boxes[1:])):
        return False
    for label in labels:
        label.set_rotation(angle)
        label.set_horizontalalignment("right")
        label.set_rotation_mode("anchor")
    return True


def plot_prediction_scatter(
    summary_rows: Sequence[Mapping[str, object]],
    deterministic_rows: Sequence[Mapping[str, object]],
    cact_labels: Sequence[str],
    run_dir: Path,
    overwrite: bool = False,
    *,
    cact_max: float | None = None,
) -> None:
    # Keep the standard-actuation plot and its weak-control zoom in separate
    # PNG files, both with linear axes.
    all_summary_rows = summary_rows
    all_deterministic_rows = deterministic_rows
    if cact_max is not None:
        summary_rows = [row for row in summary_rows if float(row["cact"]) <= cact_max]
        deterministic_rows = [
            row for row in deterministic_rows if float(row["cact"]) <= cact_max
        ]
        if not summary_rows:
            return
    else:
        # Keep the full-range figure focused on the standard
        # Cact=0.1--0.9 experiments; supplemental weak controls belong in the
        # dedicated linear-axis zoom figure.
        summary_rows = [row for row in summary_rows if float(row["cact"]) >= 0.1]
        deterministic_rows = [
            row for row in deterministic_rows if float(row["cact"]) >= 0.1
        ]
    plt = plotting_module()
    fig, ax = new_plot_canvas(plt)
    cmap = plt.get_cmap("viridis")
    displayed_cact_labels = [
        cact
        for cact in cact_labels
        if any(row["cact_label"] == cact for row in summary_rows)
    ]
    colors = cmap(np.linspace(0.05, 0.95, len(displayed_cact_labels)))
    for color, cact in zip(colors, displayed_cact_labels):
        rows = [row for row in summary_rows if row["cact_label"] == cact]
        if not rows:
            continue
        edgecolors = [
            "red" if str(row["cpert_label"]) == "0.3" else "black"
            for row in rows
        ]
        linewidths = [
            2.0 if str(row["cpert_label"]) == "0.3" else 0.7
            for row in rows
        ]
        ax.scatter(
            [row["delta_j_pred_mm"] for row in rows],
            [row["delta_j_true_mm"] for row in rows],
            s=[
                25.0 + 140.0 * (float(row["cpert"]) / 0.5) ** 0.7
                for row in rows
            ],
            color=color,
            edgecolors=edgecolors,
            linewidths=linewidths,
            alpha=0.85,
            label=rf"$C_{{\mathrm{{act}}}}={cact}$",
        )
    if deterministic_rows:
        ax.scatter(
            [row["delta_j_pred_mm"] for row in deterministic_rows],
            [row["delta_j_true_mm"] for row in deterministic_rows],
            marker="x",
            s=90,
            linewidth=2.0,
            color="black",
            label="deterministic",
            zorder=5,
        )

    all_values = [
        float(row[key])
        for row in summary_rows
        for key in ("delta_j_pred_mm", "delta_j_true_mm")
    ]
    if cact_max is not None:
        all_values.extend(
            float(row[key])
            for row in deterministic_rows
            for key in ("delta_j_pred_mm", "delta_j_true_mm")
        )
    lower, upper = min(all_values), max(all_values)
    span = upper - lower
    pad = 0.06 * (span if span > 0 else max(abs(lower), 1.0e-12))
    lower, upper = lower - pad, upper + pad
    ax.plot([lower, upper], [lower, upper], "k--", linewidth=1.2, label=r"$y=x$")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Predicted objective change\n" + r"$\Delta J_{\mathrm{pred}}$ (mm)")
    ax.set_ylabel("Nonlinear objective change\n" + r"$\Delta J_{\mathrm{true}}$ (mm)")
    if cact_max is not None:
        from matplotlib.ticker import ScalarFormatter

        for axis in (ax.xaxis, ax.yaxis):
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-3, 4))
            axis.set_major_formatter(formatter)
        rotate_overlapping_xticklabels(fig, ax)
    ax.grid(True, alpha=0.25)
    add_figure_legend(fig, ax)
    stem = "deltaJ_pred_vs_true"
    if cact_max is not None:
        stem += f"_zoom_Cact_le_{cact_max:g}"
        ax.set_title(rf"$C_{{\mathrm{{act}}}} \leq {cact_max:g}$")
    save_figure(fig, run_dir, stem, overwrite=overwrite)
    plt.close(fig)
    if cact_max is None:
        for zoom_max in (0.1, 0.001):
            plot_prediction_scatter(
                all_summary_rows, all_deterministic_rows, cact_labels, run_dir,
                overwrite=overwrite, cact_max=zoom_max,
            )


def plot_error_vs_cpert(
    summary_rows: Sequence[Mapping[str, object]],
    cpert_labels: Sequence[str],
    cact_labels: Sequence[str],
    run_dir: Path,
    overwrite: bool = False,
) -> None:
    plt = plotting_module()
    fig, ax = new_plot_canvas(plt)
    positions = np.arange(len(cpert_labels))
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(cact_labels)))
    for color, cact in zip(colors, cact_labels):
        lookup = {
            str(row["cpert_label"]): 100.0 * float(row["symmetric_normalized_error"])
            for row in summary_rows
            if row["cact_label"] == cact
        }
        ax.plot(
            positions,
            [lookup[label] for label in cpert_labels],
            marker="o",
            markersize=7.0,
            linewidth=2.2,
            color=color,
            label=rf"$C_{{\mathrm{{act}}}}={cact}$",
        )
    ax.set_xticks(positions)
    ax.set_xticklabels(
        cpert_labels, rotation=45, ha="right", rotation_mode="anchor"
    )
    ax.set_xlabel(r"Perturbation magnitude $C_{\mathrm{pert}}$")
    ax.set_ylabel("Symmetric normalized\nprediction error (%)")
    ax.grid(True, alpha=0.25)
    add_figure_legend(fig, ax)
    save_figure(fig, run_dir, "prediction_error_vs_Cpert", overwrite=overwrite)
    plt.close(fig)


def plot_error_vs_cact(
    summary_rows: Sequence[Mapping[str, object]],
    cpert_labels: Sequence[str],
    cact_labels: Sequence[str],
    run_dir: Path,
    overwrite: bool = False,
) -> None:
    """Plot finite-amplitude prediction error, with one curve per Cpert."""
    plt = plotting_module()
    fig, ax = new_plot_canvas(plt)
    positions = np.arange(len(cact_labels))
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(cpert_labels)))
    for color, cpert in zip(colors, cpert_labels):
        lookup = {
            str(row["cact_label"]): 100.0 * float(row["symmetric_normalized_error"])
            for row in summary_rows
            if row["cpert_label"] == cpert
        }
        highlighted = cpert == "0.3"
        ax.plot(
            positions,
            [lookup[label] for label in cact_labels],
            marker="o",
            markersize=9.0 if highlighted else 6.5,
            linewidth=3.2 if highlighted else 2.0,
            color="red" if highlighted else color,
            label=rf"$C_{{\mathrm{{pert}}}}={cpert}$",
            zorder=5 if highlighted else 2,
        )
    ax.set_xticks(positions)
    compact_cact_labels = []
    for label in cact_labels:
        value = float(label)
        if value == 0.00001:
            compact_cact_labels.append(r"$10^{-5}$")
        elif value == 0.00005:
            compact_cact_labels.append(r"$5\!\times\!10^{-5}$")
        elif value == 0.0001:
            compact_cact_labels.append(r"$10^{-4}$")
        elif value == 0.0005:
            compact_cact_labels.append(r"$5\!\times\!10^{-4}$")
        else:
            compact_cact_labels.append(label)
    ax.set_xticklabels(
        compact_cact_labels, rotation=60, ha="right", rotation_mode="anchor"
    )
    ax.set_xlabel(r"Actuation magnitude $C_{\mathrm{act}}$")
    ax.set_ylabel("Symmetric normalized\nprediction error (%)")
    ax.grid(True, alpha=0.25)
    add_figure_legend(fig, ax)
    save_figure(fig, run_dir, "prediction_error_vs_Cact", overwrite=overwrite)
    plt.close(fig)


def plot_deterministic(
    deterministic_rows: Sequence[Mapping[str, object]],
    run_dir: Path,
    overwrite: bool = False,
) -> None:
    if not deterministic_rows:
        return
    plt = plotting_module()
    cact = [row["cact"] for row in deterministic_rows]
    predicted = [row["delta_j_pred_mm"] for row in deterministic_rows]
    true = [row["delta_j_true_mm"] for row in deterministic_rows]
    fig, ax = new_plot_canvas(plt)
    ax.plot(cact, predicted, "o-", label=r"$\Delta J_{\mathrm{pred}}$")
    ax.plot(cact, true, "s-", label=r"$\Delta J_{\mathrm{true}}$")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel(r"Actuation magnitude $C_{\mathrm{act}}$")
    ax.set_ylabel("Deterministic objective\nchange (mm)")
    ax.grid(True, alpha=0.25)
    add_figure_legend(fig, ax)
    save_figure(fig, run_dir, "deterministic_deltaJ_vs_Cact", overwrite=overwrite)
    plt.close(fig)


def configuration_dict(
    args: argparse.Namespace,
    run_dir: Path,
    cache_dir: Path,
    cache_status: str,
) -> Dict[str, object]:
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_directory": str(run_dir),
        "ad_root": str(args.ad_root.resolve()),
        "nl_root": str(args.nl_root.resolve()),
        "cache_directory": str(cache_dir),
        "cache_status": cache_status,
        "workers": args.workers,
        "reference_deltap": str(args.reference_deltap.resolve()),
        "cpert": list(args.cpert),
        "cact": list(args.cact),
        "members": args.members,
        "target_row_zero_based": args.target_row,
        "target_col_zero_based": args.target_col,
        "rain_time_index_zero_based": args.rain_time_index,
        "adjoint_time_index": args.adjoint_time_index,
        "raw_to_objective_scale": args.raw_to_objective_scale,
        "error_floor_mm": args.error_floor_mm,
        "zero_tolerance": args.zero_tolerance,
        "objective_kernel_normalized": KERNEL.tolist(),
        "primary_prediction_definition": (
            "mean_k[(-1/60) * mean_l(A_QVAPOR_l) dot delta_x_actual_k]"
        ),
        "member_specific_prediction_retained_as_diagnostic": True,
        "input_files_opened_read_only": True,
        "overwrite_existing_results": args.overwrite,
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.ad_root = args.ad_root.expanduser().resolve()
    args.nl_root = args.nl_root.expanduser().resolve()
    args.output_root = args.output_root.expanduser().resolve()
    args.cache_root = args.cache_root.expanduser().resolve()
    args.reference_deltap = args.reference_deltap.expanduser().resolve()

    print(f"AD root (read-only): {args.ad_root}")
    print(f"NL root (read-only): {args.nl_root}")
    if args.dry_run:
        print("Checking complete paired input inventory...", flush=True)
        validate_input_paths(args)
        print("Input inventory is complete.")
        print("Dry run complete; no files or directories were created.")
        return 0

    ensure_output_is_separate(args.output_root, (args.ad_root, args.nl_root))
    requested_run_dir = args.output_root / (args.run_name or DEFAULT_RUN_NAME)
    if requested_run_dir.exists():
        if not requested_run_dir.is_dir():
            raise FileExistsError(
                f"Result path exists and is not a directory: {requested_run_dir}"
            )
        if not args.overwrite:
            raise FileExistsError(
                f"Result directory already exists: {requested_run_dir}; pass --overwrite to replace its outputs"
            )

    load_scientific_dependencies()
    ensure_output_is_separate(args.cache_root, (args.ad_root, args.nl_root))
    cache, cache_dir, cache_status = obtain_cache(args)
    member_rows, summary_rows, deterministic_rows = analyze_from_cache(args, cache)
    reference_report = add_reference_deltap(summary_rows, args.reference_deltap)
    ensemble_fit = linear_fit_metrics(summary_rows)
    deterministic_fit = linear_fit_metrics(deterministic_rows)

    ensure_output_is_separate(args.output_root, (args.ad_root, args.nl_root))
    try:
        run_dir, replacing_existing = new_run_directory(
            args.output_root, args.run_name, args.overwrite
        )
    except PermissionError:
        print(
            f"Permission denied while creating output under {args.output_root}; no input data "
            "were modified.",
            file=sys.stderr,
        )
        raise

    if replacing_existing:
        print(f"Overwriting existing result files in: {run_dir}")
    write_csv(run_dir / "memberwise_results.csv", member_rows, overwrite=args.overwrite)
    write_csv(run_dir / "ensemble_summary.csv", summary_rows, overwrite=args.overwrite)
    write_csv(run_dir / "deterministic_summary.csv", deterministic_rows, overwrite=args.overwrite)
    write_csv(
        run_dir / "global_metrics.csv",
        [
            {"analysis": "ensemble", **ensemble_fit},
            {"analysis": "deterministic", **deterministic_fit},
        ],
        overwrite=args.overwrite,
    )
    save_npz(
        run_dir / "adjoint_prediction_validation.npz",
        summary_rows,
        args.cpert,
        args.cact,
        overwrite=args.overwrite,
    )

    metadata = configuration_dict(args, run_dir, cache_dir, cache_status)
    metadata["reference_deltap_check"] = reference_report
    metadata["ensemble_global_fit"] = ensemble_fit
    metadata["deterministic_global_fit"] = deterministic_fit
    metadata_mode = "w" if args.overwrite else "x"
    with (run_dir / "metadata.json").open(metadata_mode, encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    if not args.no_plots:
        plot_prediction_scatter(
            summary_rows,
            deterministic_rows,
            args.cact,
            run_dir,
            overwrite=args.overwrite,
        )
        plot_error_vs_cpert(
            summary_rows,
            args.cpert,
            args.cact,
            run_dir,
            overwrite=args.overwrite,
        )
        plot_error_vs_cact(
            summary_rows,
            args.cpert,
            args.cact,
            run_dir,
            overwrite=args.overwrite,
        )
        plot_deterministic(
            deterministic_rows, run_dir, overwrite=args.overwrite
        )

    print(f"Validation outputs written to: {run_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PermissionError as exc:
        raise SystemExit(f"Stopped because of a permission error: {exc}") from exc
    except FileExistsError as exc:
        raise SystemExit(str(exc)) from exc
