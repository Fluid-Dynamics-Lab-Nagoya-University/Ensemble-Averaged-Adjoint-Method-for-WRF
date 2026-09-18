# WRF control preparation and analysis tools

This directory contains the additional preparation, analysis and plotting tools provided for the first revision of a journal paper manuscript.

This directory provides initial-condition preparation, nonlinear WRF submission, analysis and
plotting tools. It does not distribute computed WRF outputs, extracted array
caches, analysis results or result figures. Run the model and analysis steps below
to generate these outputs from the inputs supplied in this dataset.

## Supplied inputs

Paths below are relative to the dataset root, one directory above this README.

- `clean/wrfinput_d01_clean`: unperturbed WRF initial state at 12:00 UTC on 5 July 2018.
- `pert/pert1.dat` through `pert/pert100.dat`: the fixed surface perturbation realizations
  used in these analyses; this dataset also contains additional perturbation files.
- `prepare/A_RAINNC_absG_i0M.dat`: the terminal precipitation adjoint forcing.
- `WRF_proc_AD_org/` and `WRF_proc_NL_org/`: adjoint and nonlinear run templates,
  including `wrfbdy_d01`, `namelist.input`, physical parameter tables and auxiliary files.
- `WRFPLUS-3.9.1.1_mdfv2/`: modified WRFPLUS source. Build it on the target system;
  the archive's `wrf.exe` entries are links, not a ready-to-run model executable.
- `dats/A_QVAPOR_mean_absG_i0M_NOV0.02_sg*.dat`: supplied ensemble-mean sensitivities
  used to construct the prescribed control initial conditions. They can also be
  recomputed from adjoint outputs as described below.

This dataset does not contain the full WRF output ensemble.
The dataset-root `prepare/` contains the adjoint forcing.
The template `iodir/` directories contain auxiliary/gradient files,
not a replacement for the simulation ensemble required by these analyses.

## Software

Use Python 3.10 or newer with NumPy, netCDF4 and Matplotlib. The gradient diagnostic
script uses tqdm. Install Helvetica to retain the figure typography. The supplied
MATLAB preparation scripts require MATLAB and, where `parfor` is used, Parallel
Computing Toolbox. Building/running WRF additionally requires an appropriate Fortran
compiler, MPI and NetCDF libraries. See the dataset-root `README.md` for the model
and adjoint HPC workflow.

## Prepare initial conditions

Run the following commands from this directory. `prepare_initial_conditions.py` reads
the supplied NetCDF initial state and text perturbations directly. It uses
both-sign/domain-wide actuation with a fixed norm and applies additive noise and
actuation before clipping. It modifies only surface
QVAPOR, except that clipping applies to all negative QVAPOR
values. It does not start WRF or extract results from simulations.

Inspect the inventory before preparing a large ensemble:

```bash
python3 prepare_initial_conditions.py AD --dry-run
python3 prepare_initial_conditions.py NL --dry-run
python3 prepare_initial_conditions.py NL --validation deterministic --dry-run
```

Generate the adjoint inputs in the dataset-root `prepare/`, compatible with the
dataset-root `workdir/P01_AD` through `P04_AD` scripts:

```bash
python3 prepare_initial_conditions.py AD
```

Generate controlled nonlinear inputs in separate directories for the two
validation modes:

```bash
python3 prepare_initial_conditions.py NL
python3 prepare_initial_conditions.py NL --validation deterministic
```

The default nonlinear output locations are `prepared/ensemble/` and
`prepared/deterministic/`. The default parameter grid contains all 10 `Cpert`
values and 17 `Cact` values. Ensemble preparation creates 100 inputs per standard
parameter pair; its `Cpert=0` standard inputs are identical, retained to match
the 100-member extraction interface of `validate_adjoint_prediction.py`. Weak-actuation `Cpert=0` and
deterministic validation use one input per pair. Thus the default counts are
1000 adjoint, 16208 ensemble-controlled and 170 deterministic-controlled inputs.
Preparing every input requires about 94 GB with the supplied initial-state file, so use `--cpert` and `--cact` to
prepare selected batches. Existing input files require explicit `--overwrite`.

Prepare all eight weak-actuation magnitudes in the weak-control workflow's input directory:

```bash
python3 prepare_initial_conditions.py NL --cact 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05 --output-dir adjoint_prediction_validation/weak_control_experiments/prepare
```

The included `NL_pre_weak_control.m` prepares the same eight weak-actuation
magnitudes. It checks pairing against actual adjoint baseline outputs first;
set `WRF_AD_ROOT` to their location and `NG_EQ_SG=1` or `0` for the two validation
modes. The Python preparation does not require those computed baseline outputs.
Use one preparation method, not both. A different sensitivity directory can be
passed to Python with `--sensitivity-dir`.

## Run the model and retain its outputs

Run the dataset-root `workdir/P01_AD` through `P04_AD` workflow using the prepared
adjoint inputs and `WRF_proc_AD_org/` template. Copy its collected output directories to
`wrfout/AD/` at the dataset root, preserving names such as
`woinput_absG_i0M_NOV0.02_sg_0.3/wrfout_d01_2018-07-05_120000_woinput_AD1`.
If recomputing the sensitivity fields as well, the dataset-root `visEnsAD.py` performs
raw-gradient averaging and text export. Select each `Cpert` and the corresponding
adjoint output directory in that script; export to a new directory rather than
overwriting the supplied `dats/`, then pass it via `--sensitivity-dir` before
preparing controls. Member gradients are not individually normalized before averaging.

The weak-control nonlinear job scripts are in
`adjoint_prediction_validation/weak_control_experiments/workdir/`. They follow the
`P01`--`P04` sequence and submit to the Fujitsu HPC scheduler with
`pjsub`. Their default batch covers all eight
weak-actuation magnitudes, with 7208 ensemble or 80 deterministic cases. Standard
controls use the dataset-root `workdir/` nonlinear workflow with the corresponding
parameter settings. Defaults
locate `WRF_proc_NL_org/` and `WRFPLUS-3.9.1.1_mdfv2/` relative to the dataset
root. Set `TEMPLATE_DIR`, `WRF_EXE`, `WRF_MODULE`, `RSCGRP`, `NODES`, `MPI_PROC`,
`ELAPSE` and `CASES_PER_JOB` for the target machine. Submission requires
`--submit`; the default only generates job scripts.

The preparation and job scripts use the same eight-value weak-actuation grid:
100 members per perturbed parameter pair and one for `Cpert=0` or deterministic
validation. For deterministic controls use
`NG_EQ_SG=0`, or run `submit_NL_ng0_Fujitsu_HPC_series.sh`, which wraps the same
submitter. Keep ensemble and deterministic projects separate on the HPC so their
prepared inputs, logs and outputs do not overwrite each other. Each project expects
its inputs in its own `prepare/`. `P04` collects finished results under that
project's `outputdir/`.

Retain standard nonlinear outputs at dataset-root `wrfout/NL/`, preserving
`absG_i0M_NOV0.02_ng<Cpert>_ig<Cact>_sg<Cpert>/wrfout_d01_2018-07-05_120000_woinput_NL<member>`.
Deterministic directories use `ng0` and member 1. Retain weak ensemble and
deterministic outputs in separate locations and specify them with `--weak-root`
in the analysis commands below. Every output must contain hourly records from
12:00 through 18:00 UTC, not merely a WRF success message.

## Analyze WRF outputs

The default input roots are dataset-root `wrfout/AD/` and `wrfout/NL/`.
`WRF_AD_ROOT` and `WRF_NL_ROOT` can override them, including for the two weak-control
analysis scripts. The main script also accepts `--ad-root` and `--nl-root`.
For unperturbed validation, `analyze_unperturbed_weak_control.py --old-nl-root`
selects the standard-actuation nonlinear output directory, while `--weak-root`
selects the weak-actuation output directory.

```bash
python3 adjoint_prediction_validation/validate_adjoint_prediction.py
python3 adjoint_prediction_validation/analyze_with_weak_control.py --weak-root /path/to/ensemble/weak/outputdir
python3 adjoint_prediction_validation/analyze_unperturbed_weak_control.py --weak-root /path/to/deterministic/weak/outputdir
python3 ensemble_sampling_uncertainty/analyze_ensemble_subsample_50.py
python3 curve_chart_misc.py
```

The prediction scripts read actual NetCDF outputs and compute the objective and
post-clipping initial-state increments. Intermediate numerical caches are created
locally by these calculations, not distributed as inputs. The sampling script
uses the standard 100-member cache generated by `validate_adjoint_prediction.py` and chooses 5000
subsets of 50 distinct members with seed `20260904`. Different subsets may overlap;
each subset uses paired gradient, baseline and controlled member indices across
parameters. These are not independent new-seed WRF runs.

The objective is terminal accumulated precipitation in mm at zero-based target
indices `(34,16)`, with the normalized kernel

```text
[[0.0625, 0.125, 0.0625],
 [0.125,  0.25,  0.125 ],
 [0.0625, 0.125, 0.0625]]
```

The nonlinear change is controlled minus baseline precipitation. The adjoint
prediction is `(-1/60) * mean(A_QVAPOR)` dotted with the actual initial increment,
then averaged across members. Positive reduction is the negative of that change.

For the diagnostic plots, `statisticPlotMisc.py` directly reads the adjoint-gradient
NetCDF files under dataset-root `wrfout/diagnostics/`; `WRF_DIAGNOSTIC_ROOT` can
override this location. Use the AD directory/member naming convention above.
To read the same adjoint outputs used for prediction validation without copying
them, set `WRF_DIAGNOSTIC_ROOT` to the dataset-root `wrfout/AD/` directory.
The gradient metrics are computed from the NetCDF inputs. The control-performance
overlay uses fixed reference values stored in the script's `best_ctrl` dictionary;
it is not recalculated from nonlinear outputs by this script.
Run it from this directory; figures are written to `figs/`:

```bash
MPLBACKEND=Agg python3 statisticPlotMisc.py
```

Result images and tables are generated by running the tools; they are not
included in this supplementary directory. Calculation and submission scripts do
not install dependencies or submit jobs without the corresponding command. See
the dataset-root `LICENSE.txt` and bundled third-party licenses for model/source licensing.
