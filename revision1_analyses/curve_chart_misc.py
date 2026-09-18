#!/usr/bin/env python3
"""Plot control-response curves from nonlinear objective-change summaries.

by Shan Jiang, FDL, Nagoya University
"""

import argparse
import csv
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'adjoint_prediction_validation' / 'results'
UNPERTURBED = RESULTS / 'unperturbed_with_weak_control' / 'unperturbed_summary.csv'
ENSEMBLE = RESULTS / 'main' / 'ensemble_summary.csv'


def read_curves(path, expected_members):
    """Read finite-amplitude NL reductions; omit every supplemental Cact < 0.1."""
    curves = {}
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        required = {'cpert', 'cact', 'n_members', 'delta_j_true_mm'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: missing columns {sorted(required - set(reader.fieldnames or []))}')
        for row in reader:
            cp, ca = float(row['cpert']), float(row['cact'])
            if not math.isfinite(cp) or not math.isfinite(ca):
                raise ValueError(f'{path}: non-finite parameter')
            if ca < 0.1:
                continue
            reduction = -float(row['delta_j_true_mm'])
            if not math.isfinite(reduction):
                raise ValueError(f'{path}: non-finite response at Cact={ca}')
            if int(row['n_members']) != expected_members:
                raise ValueError(f'{path}: unexpected validation member count at Cact={ca}')
            curve = curves.setdefault(cp, {})
            if ca in curve:
                raise ValueError(f'{path}: duplicate Cpert={cp}, Cact={ca}')
            curve[ca] = reduction
    if not curves:
        raise ValueError(f'{path}: no Cact>=0.1 results')
    return curves


def load_curves():
    deterministic = read_curves(UNPERTURBED, 1)
    ensemble = read_curves(ENSEMBLE, 100)
    if set(deterministic) != set(ensemble):
        raise ValueError('The two validation summaries have different Cpert grids.')
    all_curves = list(deterministic.values()) + list(ensemble.values())
    if any(set(curve) != set(all_curves[0]) for curve in all_curves[1:]):
        raise ValueError('The curves have different Cact grids; regenerate the result summaries.')
    return sorted(all_curves[0]), deterministic, ensemble


def self_check():
    """Small fixture check for Cpert grouping, filtering, sign and validation."""
    import tempfile

    with tempfile.TemporaryDirectory(prefix='figure5-check-') as directory:
        path = Path(directory) / 'summary.csv'
        path.write_text(
            'cpert,cact,n_members,delta_j_true_mm\n'
            '0.3,0.00001,100,-999\n'
            '0,0.1,100,-999\n'
            '0.3,0.9,100,-3\n'
            '0.3,0.1,100,-2\n', encoding='utf-8',
        )
        curves = read_curves(path, 100)
        assert curves == {0.0: {0.1: 999.0}, 0.3: {0.1: 2.0, 0.9: 3.0}}
        assert all(0.0 not in curve for curve in curves.values())
        for extra in ('0.3,0.1,100,-2\n', '0.3,0.2,100,nan\n', '0.3,0.2,1,-1\n'):
            original = path.read_text(encoding='utf-8')
            path.write_text(original + extra, encoding='utf-8')
            try:
                read_curves(path, 100)
            except ValueError:
                pass
            else:
                raise AssertionError('Invalid or duplicate input was accepted')
            path.write_text(original, encoding='utf-8')
    print('Self-check passed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'valid_curve.png')
    parser.add_argument('--overwrite', action='store_true', help='Replace the two validation figures if they already exist.')
    parser.add_argument('--self-check', action='store_true')
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    output = args.output.expanduser().resolve()
    if output.suffix.lower() not in ('.png', '.jpg', '.jpeg', '.pdf', '.svg'):
        parser.error('--output must be a PNG, JPEG, PDF or SVG figure')
    ensemble_output = output.with_name(output.stem + '_ensemble' + output.suffix)
    for destination in (output, ensemble_output):
        if destination.exists() and not args.overwrite:
            parser.error(f'Figure exists: {destination}; use --overwrite to replace it')
    x, deterministic, ensemble = load_curves()

    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica'],
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Helvetica',
        'mathtext.it': 'Helvetica:italic',
        'mathtext.bf': 'Helvetica:bold',
        'font.size': 26,
        'axes.labelsize': 28,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
    })
    cpert_values = sorted(deterministic)
    values = [value for curves in (deterministic, ensemble)
              for curve in curves.values() for value in curve.values()]
    padding = 0.05 * (max(values) - min(values)) or 0.1
    colors = plt.get_cmap('cool')
    markers = ['s', 'o', 'D', 'v', 'P', 'X']
    for curves, title, destination in (
        (deterministic, 'Deterministic validation', output),
        (ensemble, 'Ensemble validation', ensemble_output),
    ):
        fig, axis = plt.subplots(figsize=(6.4, 8.2))
        axis.set_box_aspect(1)
        for index, cp in enumerate(cpert_values):
            color = colors(0.1 + 0.8 * index / max(1, len(cpert_values) - 1))
            marker, linewidth, zorder = markers[index % len(markers)], 1.4, 2
            if cp == 0.0:
                color, marker, linewidth, zorder = 'black', '^', 2.4, 4
            elif math.isclose(cp, 0.3, rel_tol=0, abs_tol=1e-12):
                color, marker, linewidth, zorder = '#FF0000', 'o', 2.4, 5
            axis.plot(x, [curves[cp][ca] for ca in x], marker=marker, color=color,
                      label=f'{cp:g}', linewidth=linewidth, markersize=5, zorder=zorder)
        axis.set_title(title, fontsize=26, pad=14)
        axis.set_xlabel(r'$C_{\mathrm{act}}$', labelpad=8)
        axis.set_ylabel(r'$-\Delta\mathcal{J}_{\mathrm{NL}}$ (mm)', labelpad=8)
        axis.set_xticks(x)
        axis.set_xticklabels([f'{ca:g}' for ca in x])
        axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
        axis.set_ylim(min(values) - padding, max(values) + padding)
        axis.grid(True, alpha=0.6)
        axis.margins(x=0.025)
        fig.subplots_adjust(left=0.23, right=0.97, bottom=0.31, top=0.94)
        fig.legend(*axis.get_legend_handles_labels(), loc='lower center',
                   bbox_to_anchor=(0.55, 0.015), fontsize=19, frameon=False,
                   title=r'$C_{\mathrm{pert}}$', title_fontsize=24, ncol=5,
                   columnspacing=0.9, handlelength=1.5, handletextpad=0.4)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination, dpi=600, bbox_inches='tight')
        plt.close(fig)
        print(f'Saved: {destination}')
    print(f'Read deterministic validation: {UNPERTURBED}')
    print(f'Read 100-member ensemble validation: {ENSEMBLE}')
    print(f'Cpert: {", ".join(f"{cp:g}" for cp in cpert_values)}')
    print('Only existing Cact>=0.1 results are plotted; no zero-control point or extrapolation.')


if __name__ == '__main__':
    main()
