#!/bin/bash
#PJM -L rscgrp=fx-small
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=05:00:00
#PJM -j
#PJM -N deltaP_npy

set -euo pipefail
trap 'echo "[ERR] line=$LINENO cmd=$BASH_COMMAND" >&2' ERR

module load netcdf-fortran
module load netcdf-c

# Set these three paths to match your HPC directory layout.
# NL_ROOT: directory containing all absG_i0M_NOV0.02_ng*_ig*_sg* subdirectories (NL outputdir)
# REF_ROOT: directory containing all woinput_absG_i0M_NOV0.02_sg_* subdirectories (AD outputdir)
# OUT_DIR:  output directory for npy files (must already exist)
export NL_ROOT="your_NL_outputdir"
export REF_ROOT="your_AD_outputdir"
export OUT_DIR="your_npy_outputdir"

# Space-separated lists, consistent with Python; ng=sg logic is handled internally
export IG_LIST="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
export SG_LIST="0 0.001 0.01 0.05 0.1 0.15 0.2 0.3 0.4 0.5"

# Input mode switch
# SWITCH1: 1=no suffix, 2=P (positive input only), 3=N (negative input only)
# SWITCH2: 1=no actr,   2=with actr (8-point actuator)
# Example: SWITCH1=3 SWITCH2=2 -> Nactr
SWITCH1=2
SWITCH2=1

if [ "$SWITCH1" -eq 1 ]; then
  NACTR_FIELD=""
else
  case "$SWITCH1" in
    2) _prefix="P" ;;
    3) _prefix="N" ;;
    *) echo "[ERR] SWITCH1 must be 1, 2, or 3"; exit 1 ;;
  esac
  case "$SWITCH2" in
    1) NACTR_FIELD="${_prefix}" ;;
    2) NACTR_FIELD="${_prefix}actr" ;;
    *) echo "[ERR] SWITCH2 must be 1 or 2"; exit 1 ;;
  esac
fi
export NACTR_FIELD

# Center grid point (0-based, consistent with Python)
export CENTER_X=34
export CENTER_Y=16

# OpenMP settings
export OMP_NUM_THREADS=24
export OMP_PROC_BIND=true
export OMP_PLACES=cores


HERE="$(cd "$(dirname "$0")" && pwd)"
"${HERE}/deltaP_mat_from_wrfout"

echo "[DONE] deltaP npy generation finished."

