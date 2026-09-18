#!/usr/bin/env bash
# Generates and submits WRF nonlinear job scripts to the Fujitsu HPC scheduler (pjsub).
# Each job covers a batch of ensemble members for one parameter combination.
#
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
mode="${1:---generate-only}"
[[ $# -le 1 && ( "$mode" == --generate-only || "$mode" == --submit ) ]] || die 'Usage: bash P03_NL_submit_wrf_jobs_Fujitsu_HPC_series.sh [--generate-only|--submit]'
[[ "$MPI_PROC" =~ ^[1-9][0-9]*$ ]] || die 'MPI_PROC must be a positive integer'
[[ "$NODES" =~ ^[1-9][0-9]*$ ]] || die 'NODES must be a positive integer'
[[ "$CASES_PER_JOB" =~ ^[1-9][0-9]*$ ]] || die 'CASES_PER_JOB must be a positive integer'
[[ "$ELAPSE" =~ ^[0-9]+:[0-5][0-9]:[0-5][0-9]$ ]] || die 'ELAPSE must be HH:MM:SS'
[[ "$RSCGRP" =~ ^[A-Za-z0-9_-]+$ ]] || die 'Invalid resource group'
if [[ "$mode" == --submit ]]; then command -v pjsub >/dev/null || die 'pjsub unavailable'; fi
while read -r tag k; do
    d="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
    need_file "$d/wrfinput_d01"
    need_file "$d/namelist.input"
    [[ -x "$d/wrf.exe" ]] || die "Run P01b first: $d"
    [[ ! -d "$d/.weak_control_running" ]] || die "Running/locked member: $d"
done < <(all_cases)

mkdir -p "$SCRIPT_DIR/jobs"
batch=$(mktemp -d "$SCRIPT_DIR/jobs/batch_XXXXXXXX")
job_files=()
index=0
for cp in "${CPERT[@]}"; do
    n=$(member_count "$cp")
    for ca in "${CACT[@]}"; do
        tag=$(case_tag "$cp" "$ca")
        for ((start=1; start<=n; start+=CASES_PER_JOB)); do
            end=$((start + CASES_PER_JOB - 1))
            (( end > n )) && end=$n
            index=$((index + 1))
            list="$batch/members_$index.txt"
            : > "$list"
            for ((k=start; k<=end; k++)); do
                [[ -f "$SCRIPT_DIR/$tag/WRF_proc_NL$k/.weak_control_complete" ]] || printf '%s\n' "$k" >> "$list"
            done
            [[ -s "$list" ]] || continue
            job="$batch/job_$index.sh"
            {
                printf '#!/bin/bash\n'
                printf '#PJM -L rscgrp=%s\n#PJM -L node=%s\n' "$RSCGRP" "$NODES"
                printf '#PJM --mpi proc=%s\n#PJM -L elapse=%s\n' "$MPI_PROC" "$ELAPSE"
                printf '#PJM -j\n#PJM -N NL%s-%s\n\n' "$start" "$end"
                printf 'set -eo pipefail\nmodule load %q\n' "$WRF_MODULE"
                printf 'export OMP_NUM_THREADS=1\n\n'
                printf 'WORK_ROOT=%q\nCASE_TAG=%q\nCASE_LIST=%q\nMPI_PROC=%q\n\n' \
                    "$SCRIPT_DIR" "$tag" "$list" "$MPI_PROC"
                cat <<'RUNNER'
run_member() (
    k="$1"
    run_dir="$WORK_ROOT/$CASE_TAG/WRF_proc_NL$k"
    cd "$run_dir" || exit 1
    # mkdir is atomic: duplicate submissions must not run the same member.
    if ! mkdir .weak_control_running 2>/dev/null; then
        echo "Member $k locked: $run_dir" >&2; exit 1
    fi
    trap 'rmdir .weak_control_running' EXIT
    if [[ -f .weak_control_complete ]]; then
        echo "Member $k already completed; skipped."; exit 0
    fi
    # Never overwrite partial output from an earlier run.
    shopt -s nullglob
    previous=(wrfout* rsl.*)
    if [[ -e wrf.log ]]; then previous+=(wrf.log); fi
    if (( ${#previous[@]} )); then
        echo "Member $k has previous output; inspect/archive before retrying." >&2; exit 1
    fi
    set +e
    mpiexec -n "$MPI_PROC" ./wrf.exe </dev/null > wrf.log 2>&1
    rc=$?
    set -e
    output='wrfout_d01_2018-07-05_12:00:00'
    if (( rc != 0 )) || [[ ! -s "$output" ]] || \
        ! tail -n 80 rsl.error.0000 2>/dev/null | grep -q 'SUCCESS COMPLETE WRF'; then
        printf 'FAILED member=%s exit=%s\n' "$k" "$rc" >&2
        exit 1
    fi
    # Final timestamp must also be checked during analysis; success alone is not
    # sufficient if the supplied namelist was changed.
    printf 'job=%s member=%s finished=%s\n' "${PJM_JOBID:-manual}" "$k" "$(date -u +%FT%TZ)" > .weak_control_complete || exit 1
    echo "SUCCESS member=$k"
)

failed=0
while IFS= read -r member; do
    if ! run_member "$member"; then failed=1; fi
done < "$CASE_LIST"
exit "$failed"
RUNNER
            } > "$job"
            bash -n "$job"
            job_files+=("$job")
        done
    done
done
printf 'Generated %s jobs in %s; mode=%s\n' "${#job_files[@]}" "$batch" "$mode"
if [[ "$mode" == --submit ]]; then
    for job in "${job_files[@]}"; do
        # Submit from the batch directory so scheduler output is kept with this batch.
        id=$(cd "$batch" && pjsub "$job")
        printf '%s\t%s\n' "$id" "$job" | tee -a "$batch/pjsub_ids.tsv"
    done
fi
