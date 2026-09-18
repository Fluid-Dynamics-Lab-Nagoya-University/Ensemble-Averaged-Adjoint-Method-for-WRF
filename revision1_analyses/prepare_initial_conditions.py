"""Prepare WRF NetCDF initial conditions from the supplied dataset inputs.

by Shan Jiang, FDL, Nagoya University
"""

import argparse
import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CPERT = ['0', '0.001', '0.01', '0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5']
CACT = ['0.00001', '0.00005', '0.0001', '0.0005', '0.001', '0.005', '0.01', '0.05',
        '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']


def modified_qvapor(q, noise, amplitude, actuation=None):
    result = q.copy()
    # Match the two sequential assignments to the original QVAPOR data type.
    surface = (q[0, 0].astype(np.float64) + noise * amplitude * 0.02).astype(q.dtype)
    if actuation is not None:
        surface = (surface.astype(np.float64) + actuation).astype(q.dtype)
    result[0, 0] = surface
    return np.maximum(result, 0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['AD', 'NL'])
    p.add_argument('--validation', choices=['ensemble', 'deterministic'], default='ensemble')
    p.add_argument('--cpert', nargs='+', choices=CPERT, default=CPERT)
    p.add_argument('--cact', nargs='+', choices=CACT, default=CACT)
    p.add_argument('--members', type=int, default=100)
    p.add_argument('--sensitivity-dir', type=Path, default=ROOT / 'dats')
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--self-check', action='store_true')
    args = p.parse_args()
    if args.self_check:
        q = np.full((1, 2, 2, 2), 0.01, dtype=np.float32)
        noise = np.array([[-1, 1], [0, 0]], dtype=np.float64)
        baseline = modified_qvapor(q, noise, 1)
        controlled = modified_qvapor(q, noise, 1, np.full((2, 2), 0.02))
        assert baseline[0, 0, 0, 0] == 0
        assert np.isclose(controlled[0, 0, 0, 0], 0.01)
        assert np.array_equal(controlled[0, 1], q[0, 1])
        assert controlled.dtype == q.dtype
        print('Self-check passed: additive noise, control-before-clipping and unchanged upper layers.')
        return
    if not 1 <= args.members <= 100:
        p.error('--members must be between 1 and 100')
    if len(set(args.cpert)) != len(args.cpert) or len(set(args.cact)) != len(args.cact):
        p.error('Parameter lists must not contain duplicates')
    if args.mode == 'AD' and args.validation != 'ensemble':
        p.error('--validation applies only to NL preparation')
    default = ROOT / 'prepare' if args.mode == 'AD' else HERE / 'prepared' / args.validation
    output = (args.output_dir or default).expanduser().resolve()
    clean = ROOT / 'clean' / 'wrfinput_d01_clean'
    with Dataset(clean) as ds:
        q = np.asarray(ds['QVAPOR'][:])
    if q.shape != (1, 31, 50, 50):
        raise ValueError(f'Unexpected QVAPOR dimensions: {q.shape}')
    noises = []
    for member in range(1, args.members + 1):
        # MATLAB reshapes column-major [x,y]; NetCDF Python indexes [y,x].
        noise = np.loadtxt(ROOT / 'pert' / f'pert{member}.dat')
        if noise.size != 2500 or not np.isfinite(noise).all():
            raise ValueError(f'Invalid perturbation member {member}')
        noises.append(noise.reshape(50, 50))
    tasks, forcings = [], []
    for cp in args.cpert:
        if args.mode == 'AD':
            folder = output / f'woinput_absG_i0M_NOV0.02_sg_{cp}'
            forcings.append(folder / 'A_RAINNC_absG_i0M.dat')
            for member in range(1, args.members + 1):
                tasks.append((folder / f'wrfinput_d01_woinput_{member}', cp, member, None))
            continue
        raw = np.loadtxt(args.sensitivity_dir / f'A_QVAPOR_mean_absG_i0M_NOV0.02_sg{cp}.dat')
        if raw.size != 31 * 50 * 50 or not np.isfinite(raw).all():
            raise ValueError(f'Invalid sensitivity for Cpert={cp}')
        gradient = raw.reshape((31, 50, 50), order='F')[0]
        maximum = np.max(np.abs(gradient))
        if maximum == 0:
            raise ValueError(f'Zero sensitivity for Cpert={cp}')
        direction = gradient / maximum
        # Both-sign, domain-wide fixed-norm actuation.
        direction = direction * 0.02
        direction = direction * (0.095222 / np.linalg.norm(direction))
        ng = cp if args.validation == 'ensemble' else '0'
        for ca in args.cact:
            folder = output / f'absG_i0M_NOV0.02_ng{ng}_ig{ca}_sg{cp}'
            members = args.members
            if args.validation == 'deterministic' or (cp == '0' and float(ca) < 0.1):
                members = 1
            for member in range(1, members + 1):
                tasks.append((folder / f'wrfinput_d01_woinput_ng_{member}', ng, member,
                              direction * float(ca)))
    for destination in [t[0] for t in tasks] + forcings:
        if destination.exists() and not args.overwrite:
            raise FileExistsError(f'{destination} exists; use --overwrite')
    if args.mode == 'AD' and not (ROOT / 'prepare' / 'A_RAINNC_absG_i0M.dat').is_file():
        raise FileNotFoundError('Missing dataset-root prepare/A_RAINNC_absG_i0M.dat')
    print(f'{len(tasks)} initial-condition files; output: {output}')
    print('Estimated input storage: %.1f GB' % (len(tasks) * clean.stat().st_size / 1e9))
    if args.dry_run:
        print('Dry run: no files written and no WRF simulations launched.')
        return
    for destination in forcings:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / 'prepare' / 'A_RAINNC_absG_i0M.dat', destination)
    for count, (destination, ng, member, actuation) in enumerate(tasks, 1):
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(clean, destination)
        with Dataset(destination, 'r+') as ds:
            ds['QVAPOR'][:] = modified_qvapor(q, noises[member - 1], float(ng), actuation)
        if count % 100 == 0 or count == len(tasks):
            print(f'Prepared {count}/{len(tasks)}', flush=True)


if __name__ == '__main__':
    main()
