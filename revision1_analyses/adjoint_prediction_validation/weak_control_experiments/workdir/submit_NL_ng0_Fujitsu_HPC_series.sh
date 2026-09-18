#!/usr/bin/env bash
# Ensemble-mean sensitivity applied to the unperturbed initial state.
# Reuse the P03 submitter without changing WRF settings.
#
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
export NG_EQ_SG=0
exec bash "$(dirname "$0")/P03_NL_submit_wrf_jobs_Fujitsu_HPC_series.sh" "$@"
