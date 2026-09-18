#!/usr/bin/env bash
# Shared by P01--P04. Paths may be overridden in the environment.
#
# by Shan Jiang, FDL, Nagoya University

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATASET_ROOT="$(cd "${PROJECT_DIR}/../../.." && pwd)"
TEMPLATE_DIR="${TEMPLATE_DIR:-${DATASET_ROOT}/WRF_proc_NL_org}"
WRF_EXE="${WRF_EXE:-${DATASET_ROOT}/WRFPLUS-3.9.1.1_mdfv2/main/wrf.exe}"
WRF_MODULE="${WRF_MODULE:-wrf/3.9.1.1}"
RSCGRP="${RSCGRP:-fx-small}"
NODES="${NODES:-3}"
MPI_PROC="${MPI_PROC:-96}"
ELAPSE="${ELAPSE:-48:00:00}"
CASES_PER_JOB="${CASES_PER_JOB:-20}"
CPERT=(0 0.001 0.01 0.05 0.1 0.15 0.2 0.3 0.4 0.5)
CACT=(0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05)
# Matches MATLAB ngEQsg: 1 = perturbed ensemble, 0 = original initial state.
NG_EQ_SG="${NG_EQ_SG:-1}"
if [[ "$NG_EQ_SG" != 0 && "$NG_EQ_SG" != 1 ]]; then
    echo 'ERROR: NG_EQ_SG must be 0 or 1' >&2
    exit 1
fi
case_tag() {
    local ng="$1"
    if [[ "$NG_EQ_SG" == 0 ]]; then ng=0; fi
    printf 'absG_i0M_NOV0.02_ng%s_ig%s_sg%s' "$ng" "$2" "$1"
}
member_count() {
    if [[ "$NG_EQ_SG" == 0 || "$1" == 0 ]]; then echo 1; else echo 100; fi
}
die() { echo "ERROR: $*" >&2; exit 1; }
need_file() { [[ -f "$1" ]] || die "Missing file: $1"; }
# Eight weak actuations: 7208 ensemble or 80 deterministic cases.
all_cases() {
    local cp ca k n tag
    for cp in "${CPERT[@]}"; do
        n=$(member_count "$cp")
        for ca in "${CACT[@]}"; do
            tag=$(case_tag "$cp" "$ca")
            for ((k=1; k<=n; k++)); do printf '%s %d\n' "$tag" "$k"; done
        done
    done
}
