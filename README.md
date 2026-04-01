# Dataset for ensemble-averaged adjoint sensitivity analysis validated by WRF

## Workflow

```
====local====
AD_pre_g_series.m

manual: copy prepare/ and WRF_proc_AD_org/ to HPC (AD)

====HPC AD workdir====
P01_AD → P01b_AD → P02_AD → P03_AD → [WRF adjoint runs] → P04_AD

manual: copy HPC AD outputdir/ to local wrfout/

====local====
visEnsAD.py → NL_pre_ngv_ig_sgv_series.m

manual: copy prepare/ and WRF_proc_NL_org/ to HPC (NL)

====HPC NL workdir====
P01_NL → P01b_NL → P02_NL → P03_NL → [WRF nonlinear runs] → P04_NL

====HPC NL postpro/====
build and run → deltaP_mat_*.npy

manual: copy HPC NL outputdir/ to local wrfout/  (needed for visEnsNLMisc.py)
          copy postpro npy to local dats/           (needed for cont_ensAD_par.py run_mode=2)

====local====
cont_ensAD_par.py, visEnsNLMisc.py
```

> **Note:** The HPC P (P01, P02, ...) scripts (in `workdir/`) are written for the Fujitsu HPC environment (for example, the supercomputer "FLOW" of Nagoya university) and must be run from the `workdir/` directory for the respective HPC projects (NL or AD). Users running on a different system will need to adapt the WRF source code, job scheduler and MPI settings accordingly. File transfers between local and HPC (marked as `manual` above) are performed manually by the user.

## Scripts

> **Note:** The scripts are not arranged based on execution order. See "workflow" above.

### Adjoint (AD) runs

- **`AD_pre_g_series.m`** — Generates perturbed WRF initial conditions for adjoint runs by adding Gaussian noise to QVAPOR across a range of noise rates.
- **`P01_AD_cp_WRF_proc_series.sh`** — Sets up base WRF working directories on the AD HPC project and distributes the perturbed initial conditions across ensemble members.
- **`P01b_AD_relink_wrf_series.sh`** — Re-links `wrf.exe` into each ensemble member's working directory after P01_AD. Requires setting `your_wrf_directory` to the parent folder of the local WRF installation (compiled from WRFPLUS-3.9.1.1_mdfv2).
- **`P02_AD_cp_wrfinput_series.sh`** — Distributes perturbed WRF initial conditions and adjoint forcing files (A_RAINNC) into each ensemble member's working directory on HPC.
- **`P03_AD_submit_wrf_jobs_Fujitsu_HPC_series.sh`** — Generates and submits WRF adjoint job scripts to the Fujitsu HPC scheduler. Each job covers a batch of ensemble members for one noise rate. Users may edit or recreate this script for corresponding HPC systems.
- **`P04_AD_gather_result_series.sh`** — Collects WRF adjoint output files from each ensemble member's working directory and gathers them into a single output directory per noise rate.

### Nonlinear (NL) runs

- **`NL_pre_ngv_ig_sgv_series.m`** — Generates perturbed WRF initial conditions for nonlinear runs, combining adjoint sensitivity fields with background noise across a range of input energy and noise rates.
- **`NL_pre_LA_ngv_ig_sgv_series.m`** — Variant of `NL_pre_ngv_ig_sgv_series.m` that uses the linear-approximation sensitivity field (`Smat` from `linear_approx_pre.m`) instead of the adjoint field. Must be run in the same MATLAB session as `linear_approx_pre.m`. For validating pseudo-inverse results.
- **`P01_NL_cp_WRF_proc_series.sh`** — Sets up base WRF working directories on the NL HPC project and distributes the perturbed initial conditions across ensemble members.
- **`P01b_NL_relink_wrf_series.sh`** — Re-links `wrf.exe` into each ensemble member's working directory after P01_NL. Requires setting `your_wrf_directory` to the parent folder of the local WRF installation (compiled from WRFPLUS-3.9.1.1_mdfv2).
- **`P02_NL_cp_wrfinput_series.sh`** — Distributes perturbed WRF initial conditions into each ensemble member's working directory on HPC.
- **`P03_NL_submit_wrf_jobs_Fujitsu_HPC_series.sh`** — Generates and submits WRF nonlinear job scripts to the Fujitsu HPC scheduler. Each job covers a batch of ensemble members for one parameter combination.
- **`P04_NL_gather_result_series.sh`** — Collects WRF nonlinear output files from each ensemble member's working directory and gathers them into a single output directory per parameter combination.

### Post-processing (HPC — NL project, postpro/)

- **`deltaP_mat_from_wrfout.f90`** — Fortran program that reads NL and AD WRF output on the HPC, computes ensemble-averaged ΔP across all parameter combinations, and writes `.npy`-compatible output.
- **`build_deltaP_npy.sh`** — Compiles and runs `deltaP_mat_from_wrfout.f90` locally on the HPC node.
- **`submit_deltaP_npy.pjm.sh`** — Submits the Fortran post-processing job to the Fujitsu HPC scheduler.

### Analysis and visualization (local)

- **`linear_approx_pre.m`** — Estimates the sensitivity field through linear approximation (pseudo-inverse). Constructs a perturbation matrix D from ensemble noise files, computes the precipitation response vector r using a 3×3 approximate Gaussian-weighted kernel at a target grid point, and recovers the sensitivity as s = r′ · pinv(D). Produces diagnostic figures of the noise pattern, clean RAINNC, clean QVAPOR, and the estimated sensitivity field. Should be used together with `NL_pre_LA_ngv_ig_sgv_series.m`.
- **`visEnsAD.py`** — Reads WRF adjoint output, computes the ensemble average of the sensitivity field (A_QVAPOR), saves the result in 2500*1 (depends on simulation resolution) vectors as a `.dat` file, and produces a map figure.
- **`cont_ensAD_par.py`** — Reads NL and AD WRF output (or pre-computed `.npy`), computes the ensemble-averaged precipitation change (ΔP) across a grid of input energy and noise rates, saves a ΔP matrix as `.npy` (if not pre-computed), and plots a contour figure.
- **`visEnsNLMisc.py`** — Produces two types of diagnostic figures from NL output: (mode 1) QVAPOR change ratio map derived from the ensemble-averaged sensitivity field; (mode 2) ensemble-averaged RAINNC difference between NL and AD runs.
- **`curve_chart.py`** — Reads ΔP matrices from `.npy` files and plots precipitation reduction curves as a function of input energy rate for different sensitivity/validation configurations.
- **`schematic_landscape_adv.py`** — Toy model illustrating the ensemble-averaged adjoint concept. Generates schematic figures of a synthetic cost-function landscape, showing how ensemble-averaged sensitivity directions compare to deterministic and ideal directions across varying ensemble sizes and perturbation magnitudes. This file is independent of the WRF simulations.
- **`schematic_landscape_3D.py`** — Another toy-model script using the same synthetic landscape as `schematic_landscape_adv.py`. Produces a 3D surface plot of the cost-function landscape with a dashed arrow from the starting point to the global minimum. Independent of the WRF simulations.
- **`statisticPlotMisc.py`** — Reads WRF adjoint output across all noise rates, computes coherence and deflection from the ensemble gradient fields, and produces two figures: (1) deflection and coherence trends versus C_pert; (2) normalized composite metrics overlaid with the control performance curve.
- **`func_analysis.py`** — Helper
- **`plotting_tools.py`** — Helper

## Source code

- **`WRFPLUS-3.9.1.1_mdfv2/`** — Modified WRFPLUS 3.9.1.1 source code enabling standalone adjoint (AD) computation.

## File structure

### Local

```
.
├── WRFPLUS-3.9.1.1_mdfv2/          # Modified WRFPLUS with standalone adjoint (AD) runs.
├── AD_pre_g_series.m
├── linear_approx_pre.m
├── NL_pre_LA_ngv_ig_sgv_series.m
├── visEnsAD.py
├── NL_pre_ngv_ig_sgv_series.m
├── cont_ensAD_par.py
├── visEnsNLMisc.py
├── curve_chart.py
├── schematic_landscape_adv.py
├── schematic_landscape_3D.py
├── statisticPlotMisc.py
├── func_analysis.py
├── plotting_tools.py
├── customcmap.rgb
├── workdir/                        # HPC P scripts (AD and NL); copy to HPC workdir/ before running
├── postpro/                        # HPC-side Fortran post-processing; copy to HPC NL project
├── WRF_proc_AD_org/                # WRF adjoint template directory (copy to HPC manually)
├── WRF_proc_NL_org/                # WRF nonlinear template directory (copy to HPC manually)
├── clean/                          # base WRF input file without noise and actuation
├── pert/                           # pre-generated perturbation patterns
├── prepare/                        # perturbed initial conditions (copy to HPC manually)
├── wrfout/                         # WRF output (AD and NL; manually copy from HPC)
├── dats/                           # data for analysis
└── figs/                           # figures
```

### HPC — AD project directory

```
.
├── WRF_proc_AD_org/                # WRF template directory (copy from local)
├── prepare/                        # copy from local
├── inputdir/                       # generated by P01_AD
├── workdir/                        # all AD P scripts run from here; generated by P01_AD
└── outputdir/                      # WRF adjoint output collected by P04_AD
```

### HPC — NL project directory

```
.
├── WRF_proc_NL_org/                # WRF template directory (copy from local)
├── prepare/                        # copy from local
├── inputdir/                       # generated by P01_NL
├── workdir/                        # all NL P scripts run from here; generated by P01_NL
├── outputdir/                      # WRF nonlinear output collected by P04_NL
└── postpro/                        # copy from local; produces deltaP_mat_*.npy
```
