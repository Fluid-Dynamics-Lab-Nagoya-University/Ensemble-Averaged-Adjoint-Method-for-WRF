#!/usr/bin/env python3
"""Compare adjoint predictions and nonlinear responses at 17 actuation magnitudes.

Reads standard and weak-control WRF outputs without modifying them or running WRF.
Replacing existing result files requires --overwrite.

by Shan Jiang, FDL, Nagoya University
"""
from __future__ import annotations
import argparse
import copy
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import validate_adjoint_prediction as base

HERE = Path(__file__).resolve().parent
# Supplemental weak-actuation simulations stored under weak_control_experiments.
# Explicit decimal labels match the shell and MATLAB directory names.
ULTRAWEAK = [
    '0.00001',
    '0.00005',
    '0.0001',
    '0.0005',
    '0.001',
]
WEAK = ULTRAWEAK + [
    '0.005',
    '0.01',
    '0.05',
]


def check_file(path_string):
    """Check Times and variable dimensions; no full simulation arrays are read."""
    base.load_scientific_dependencies()
    with base.Dataset(path_string, 'r') as ds:
        for name in ['Times', 'QVAPOR', 'RAINNC']:
            if name not in ds.variables:
                raise ValueError(f'{path_string}: missing {name}')
        q, rain = ds['QVAPOR'], ds['RAINNC']
        if q.ndim != 4 or rain.ndim != 3 or q.shape[1] < 1:
            raise ValueError(f'{path_string}: invalid variable dimensions')
        if q.shape[-2:] != (50, 50) or rain.shape[-2:] != (50, 50):
            raise ValueError(f'{path_string}: expected 50x50 grid')
        times = [b''.join(row).decode('ascii').strip('\x00 ') for row in ds['Times'][:]]
        expected = [f'2018-07-05_{hour:02d}:00:00' for hour in range(12, 19)]
        if times != expected or q.shape[0] != 7 or rain.shape[0] != 7:
            raise ValueError(f'{path_string}: expected hourly 12:00--18:00, got {times}')
    return path_string


def select_old(cache, start, stop, members):
    result = {}
    for key, values in cache.items():
        if key in ('s_raw', 'j_baseline_mm'):
            result[key] = values[start:stop, :members]
        else:
            result[key] = values[start:stop, :, :members]
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--weak-root', type=Path, default=HERE/'weak_control_experiments/outputdir')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--run-name', default='with_weak_control')
    p.add_argument('--overwrite', action='store_true', help='Allow existing outputs in the selected result directory to be replaced.')
    p.add_argument('--check-only', action='store_true', help='Validate weak-control NetCDF files and Times, write nothing.')
    p.add_argument('--no-plots', action='store_true')
    p.add_argument('--from-cache', action='store_true', help='Require existing caches; skip all NetCDF reads and inventory checks.')
    opt = p.parse_args()
    if opt.from_cache and opt.check_only:
        p.error('--from-cache cannot be combined with --check-only')
    if opt.workers < 1:
        p.error('--workers must be positive')
    if Path(opt.run_name).name != opt.run_name or opt.run_name in ('', '.', '..'):
        p.error('--run-name must be a single directory name')
    opt.weak_root = opt.weak_root.expanduser().resolve()
    old_args = base.parse_args([])
    old_args.workers = opt.workers
    old_args.from_cache = opt.from_cache
    old_args.overwrite = opt.overwrite
    out = HERE/'results'/opt.run_name
    base.ensure_output_is_separate(out, (old_args.ad_root, old_args.nl_root, opt.weak_root))
    if not opt.check_only and out.exists() and not opt.overwrite:
        p.error(f'Result directory exists: {out}; pass --overwrite to replace its outputs')
    if opt.check_only:
        paths = [base.nl_file(opt.weak_root, cp, ca, k)
                 for cp in base.DEFAULT_CPERT for ca in WEAK
                 for k in range(1, 2 if cp == '0' else 101)]
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(f'Missing {len(missing)} files; first 10:\n'+'\n'.join(missing[:10]))
        print(f'Checking {len(paths)} weak-control NL files and hourly timestamps...', flush=True)
        with ProcessPoolExecutor(max_workers=opt.workers,
                                 mp_context=multiprocessing.get_context('spawn')) as pool:
            for count, _ in enumerate(pool.map(check_file, map(str, paths)), 1):
                if count % 100 == 0:
                    print(f'Checked {count}/{len(paths)}', flush=True)
        print('Weak-control input checks passed.', flush=True)
        print('Check only: no files written; no prediction analysis performed.')
        return
    base.load_scientific_dependencies()
    old, old_cache_path, old_cache_status = base.obtain_cache(old_args)
    member_rows, summary_rows, deterministic_rows, provenance = [], [], [], []
    # Recompute standard-actuation summaries, with only member 1 for Cpert=0.
    for start, stop, members in [(0, 1, 1), (1, 10, 100)]:
        args = copy.deepcopy(old_args)
        args.cpert = list(base.DEFAULT_CPERT[start:stop])
        args.members = members
        subset = select_old(old, start, stop, members)
        m, s, d = base.analyze_from_cache(args, subset)
        member_rows.extend(m); summary_rows.extend(s); deterministic_rows.extend(d)
        args.nl_root = opt.weak_root
        args.cact = WEAK
        args.workers = opt.workers
        args.from_cache = opt.from_cache
        args.cache_root = HERE/'cache_weak_control'
        base.ensure_output_is_separate(args.cache_root, (args.ad_root, args.nl_root))
        new, cache_path, status = base.obtain_cache(args)
        for key in ('s_raw', 'j_baseline_mm'):
            if not base.np.allclose(new[key], subset[key], rtol=0, atol=1e-8):
                raise ValueError(f'{key}: weak-control extraction differs from standard cache; check AD inputs')
        m, s, d = base.analyze_from_cache(args, new)
        member_rows.extend(m); summary_rows.extend(s); deterministic_rows.extend(d)
        provenance.append({'cache': str(cache_path), 'status': status,
                           'cpert': args.cpert, 'members': members})
    for rows in (member_rows, summary_rows):
        rows.sort(key=lambda r: (r['cpert'], r['cact'], r.get('member', 0)))
    deterministic_rows.sort(key=lambda r: r['cact'])
    for row in member_rows + summary_rows:
        pred, true = row['delta_j_pred_mm'], row['delta_j_true_mm']
        if not base.np.isfinite(pred) or not base.np.isfinite(true):
            raise ValueError('Nonfinite prediction or response encountered')
        # Add unfloored relative errors for weak responses.
        row['pred_over_true'] = pred / true if true != 0 else float('nan')
        row['truth_relative_error_without_floor'] = abs(pred-true)/abs(true) if true != 0 else float('nan')
    cact = sorted(list(base.DEFAULT_CACT)+WEAK, key=float)
    expected_summary_count = len(base.DEFAULT_CPERT) * len(cact)
    old_member_count = len(base.DEFAULT_CACT) * (1 + 100 * (len(base.DEFAULT_CPERT) - 1))
    weak_member_count = len(WEAK) * (1 + 100 * (len(base.DEFAULT_CPERT) - 1))
    expected_member_count = old_member_count + weak_member_count
    if len(summary_rows) != expected_summary_count or len(member_rows) != expected_member_count:
        raise ValueError(
            'Unexpected combined result count: '
            f'got {len(summary_rows)} summaries/{len(member_rows)} members, '
            f'expected {expected_summary_count}/{expected_member_count}'
        )
    reference = base.add_reference_deltap(summary_rows, old_args.reference_deltap)
    metrics = []
    metric_groups = [
        ('all', summary_rows),
        ('weak', [r for r in summary_rows if r['cact_label'] in WEAK]),
        ('ultraweak', [r for r in summary_rows if r['cact_label'] in ULTRAWEAK]),
    ]
    for label, rows in metric_groups:
        metrics.append({'analysis': label, **base.linear_fit_metrics(rows)})
    out.mkdir(parents=True, exist_ok=opt.overwrite)
    base.write_csv(out/'memberwise_results.csv', member_rows, overwrite=opt.overwrite)
    base.write_csv(out/'ensemble_summary.csv', summary_rows, overwrite=opt.overwrite)
    base.write_csv(out/'deterministic_summary.csv', deterministic_rows, overwrite=opt.overwrite)
    base.write_csv(out/'global_metrics.csv', metrics, overwrite=opt.overwrite)
    base.save_npz(out/'adjoint_prediction_validation.npz', summary_rows, base.DEFAULT_CPERT, cact, overwrite=opt.overwrite)
    metadata = base.configuration_dict(old_args, out, old_cache_path, old_cache_status)
    metadata.update(created_utc=datetime.now(timezone.utc).isoformat(), cact=cact,
                    members_by_cpert={cp: 1 if cp == '0' else 100 for cp in base.DEFAULT_CPERT},
                    weak_nl_root=str(opt.weak_root), weak_caches=provenance,
                    reference_deltap_check=reference,
                    note='Old Cpert=0 uses member 1; additional ratio errors have no denominator floor.')
    metadata.pop('members', None)
    metadata_mode = 'w' if opt.overwrite else 'x'
    with (out/'metadata.json').open(metadata_mode, encoding='utf-8') as handle:
        handle.write(json.dumps(metadata, indent=2)+'\n')
    if not opt.no_plots:
        base.plot_prediction_scatter(summary_rows, deterministic_rows, cact, out, overwrite=opt.overwrite)
        base.plot_error_vs_cpert(summary_rows, base.DEFAULT_CPERT, cact, out, overwrite=opt.overwrite)
        base.plot_error_vs_cact(summary_rows, base.DEFAULT_CPERT, cact, out, overwrite=opt.overwrite)
        base.plot_deterministic(deterministic_rows, out, overwrite=opt.overwrite)
    print(
        f'Completed: {len(summary_rows)} comparisons, '
        f'{len(member_rows)} member records. Results: {out}'
    )


if __name__ == '__main__':
    main()
