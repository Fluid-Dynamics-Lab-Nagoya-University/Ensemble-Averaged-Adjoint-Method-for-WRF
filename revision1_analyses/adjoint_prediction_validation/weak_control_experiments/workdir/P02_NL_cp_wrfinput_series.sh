#!/usr/bin/env bash
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
while read -r tag k; do
    need_file "$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k"
    need_file "$SCRIPT_DIR/$tag/WRF_proc_NL$k/namelist.input"
    [[ ! -e "$SCRIPT_DIR/$tag/WRF_proc_NL$k/wrfinput_d01" && ! -L "$SCRIPT_DIR/$tag/WRF_proc_NL$k/wrfinput_d01" ]] || die "Existing run input: $tag/$k"
done < <(all_cases)
while read -r tag k; do
    cp "$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k" "$SCRIPT_DIR/$tag/WRF_proc_NL$k/wrfinput_d01"
done < <(all_cases)
echo 'P02 complete.'
