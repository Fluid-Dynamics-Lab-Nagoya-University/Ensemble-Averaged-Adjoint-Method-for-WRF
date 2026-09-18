#!/usr/bin/env bash
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
need_file "$TEMPLATE_DIR/namelist.input"
# Preflight the complete inventory before copying anything.
while read -r tag k; do
    need_file "$PROJECT_DIR/prepare/$tag/wrfinput_d01_woinput_ng_$k"
    [[ ! -e "$SCRIPT_DIR/$tag/WRF_proc_NL$k" ]] || die "Existing run directory: $tag/WRF_proc_NL$k"
    [[ ! -e "$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k" ]] || die "Existing staged input: $tag/$k"
done < <(all_cases)
while read -r tag k; do
    dest="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
    mkdir -p "$dest" "$PROJECT_DIR/inputdir/$tag"
    # Exclude previous inputs/results and executable; P01b/P02 install these.
    rsync -a --exclude='wrfinput_d01' --exclude='wrf.exe' --exclude='rsl.*' \
        --exclude='wrfout*' --exclude='wrfrst*' --exclude='*.log' "$TEMPLATE_DIR/" "$dest/"
    cp "$PROJECT_DIR/prepare/$tag/wrfinput_d01_woinput_ng_$k" \
        "$PROJECT_DIR/inputdir/$tag/wrfinput_d01_woinput_$k"
done < <(all_cases)
echo "P01 complete: $(all_cases | wc -l | tr -d ' ') run directories and staged inputs."
