#!/usr/bin/env bash
# by Shan Jiang, FDL, Nagoya University

set -euo pipefail
source "$(dirname "$0")/config.sh"
[[ -x "$WRF_EXE" ]] || die "WRF executable missing/not executable: $WRF_EXE"
while read -r tag k; do
    dest="$SCRIPT_DIR/$tag/WRF_proc_NL$k"
    [[ -d "$dest" ]] || die "Run P01 first: $dest"
    [[ ! -e "$dest/wrf.exe" && ! -L "$dest/wrf.exe" ]] || die "Existing wrf.exe: $dest"
done < <(all_cases)
while read -r tag k; do
    ln -s "$WRF_EXE" "$SCRIPT_DIR/$tag/WRF_proc_NL$k/wrf.exe"
done < <(all_cases)
echo 'P01b complete.'
