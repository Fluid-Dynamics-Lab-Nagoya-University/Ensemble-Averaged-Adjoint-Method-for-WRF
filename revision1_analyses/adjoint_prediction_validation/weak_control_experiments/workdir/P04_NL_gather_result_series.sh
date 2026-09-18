#!/usr/bin/env bash
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
missing=0
collected=0
mkdir -p "$PROJECT_DIR/outputdir"
report=$(mktemp "$PROJECT_DIR/outputdir/incomplete_XXXXXXXX.txt")
while read -r tag k; do
    run_dir="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
    src="$run_dir/wrfout_d01_2018-07-05_12:00:00"
    dest="$PROJECT_DIR/outputdir/$tag/wrfout_d01_2018-07-05_120000_woinput_NL$k"
    if [[ ! -f "$run_dir/.weak_control_complete" || ! -s "$src" ]]; then
        printf '%s %s\n' "$tag" "$k" >> "$report"
        missing=$((missing + 1)); continue
    fi
    mkdir -p "$(dirname "$dest")"
    if [[ -e "$dest" ]]; then
        cmp -s "$src" "$dest" || die "Existing collected file differs: $dest"
    else
        tmp=$(mktemp "${dest}.tmp_XXXXXXXX")
        cp "$src" "$tmp"
        # Atomic, exclusive publication; never replace an existing result.
        ln "$tmp" "$dest"
        unlink "$tmp"
    fi
    collected=$((collected + 1))
done < <(all_cases)
echo "Collected/verified: $collected; incomplete: $missing; report: $report"
(( missing == 0 ))
