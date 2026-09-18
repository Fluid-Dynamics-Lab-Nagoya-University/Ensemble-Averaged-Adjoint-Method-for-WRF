#!/usr/bin/env bash
# Run directly inside an allocated interactive job; does not submit scheduler jobs.
# Same copies and exclusions as P01_NL_cp_WRF_proc_series.sh, with independent member workers.
#
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
workers="${P01_WORKERS:-16}"
[[ "$workers" =~ ^[1-9][0-9]*$ ]] || die 'P01_WORKERS must be a positive integer'
(( workers <= 16 )) || die 'This version supports at most 16 workers'
[[ $# == 0 ]] || die 'Usage: P01_WORKERS=16 bash P01_NL_cp_WRF_proc_series_parallel.sh'
need_file "$TEMPLATE_DIR/namelist.input"
command -v rsync >/dev/null || die 'rsync is unavailable'

tags=()
members=()
echo 'Checking all inputs and destinations before copying...'
while read -r tag k; do
    need_file "$PROJECT_DIR/prepare/$tag/wrfinput_d01_woinput_ng_$k"
    dest="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
    staged="$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k"
    [[ ! -e "$dest" && ! -L "$dest" ]] || die "Existing run directory: $dest"
    [[ ! -e "$staged" && ! -L "$staged" ]] || die "Existing staged input: $staged"
    tags+=("$tag")
    members+=("$k")
done < <(all_cases)
total=${#tags[@]}
(( total > 0 )) || die 'No cases configured'
(( workers <= total )) || workers=$total
mkdir -p "$SCRIPT_DIR/logs"
log_dir=$(mktemp -d "$SCRIPT_DIR/logs/P01_parallel_XXXXXXXX")
echo "Copying $total cases using $workers workers. Logs: $log_dir"

copy_lane() {
    local lane="$1" idx tag k dest staged count=0
    for ((idx=lane; idx<total; idx+=workers)); do
        tag="${tags[idx]}"; k="${members[idx]}"
        dest="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
        staged="$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k"
        echo "START $tag member=$k"
        mkdir -p "$SCRIPT_DIR/$tag" "$PROJECT_DIR/inputdir/$tag" || return 1
        # Exclusive creation also protects against two overlapping P01 invocations.
        mkdir "$dest" || return 1
        rsync -a --exclude='wrfinput_d01' --exclude='wrf.exe' --exclude='rsl.*' \
            --exclude='wrfout*' --exclude='wrfrst*' --exclude='*.log' \
            "$TEMPLATE_DIR/" "$dest/" || return 1
        cp "$PROJECT_DIR/prepare/$tag/wrfinput_d01_woinput_ng_$k" "$staged" || return 1
        count=$((count + 1))
        echo "DONE $tag member=$k"
    done
    printf '%s\n' "$count" > "$log_dir/worker_$lane.count"
}

pids=()
for ((lane=0; lane<workers; lane++)); do
    copy_lane "$lane" > "$log_dir/worker_$lane.log" 2>&1 &
    pids+=("$!")
done
failed=0
for ((lane=0; lane<workers; lane++)); do
    if ! wait "${pids[lane]}"; then
        echo "Worker $lane failed; see $log_dir/worker_$lane.log" >&2
        failed=1
    fi
done
(( failed == 0 )) || die "P01 incomplete. Inspect logs and partial directories before retrying: $log_dir"
copied=0
for ((lane=0; lane<workers; lane++)); do
    read -r count < "$log_dir/worker_$lane.count"
    copied=$((copied + count))
done
(( copied == total )) || die "Case-count mismatch: expected $total, copied $copied"
echo "P01 parallel complete: $copied run directories and staged inputs. Continue with P01b, then P02."
