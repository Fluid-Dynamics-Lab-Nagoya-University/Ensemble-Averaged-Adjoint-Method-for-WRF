#!/bin/bash
# Collects WRF adjoint output files from each ensemble member's working directory
# and gathers them into a single output directory per noise rate.
#
# by Shan Jiang, FDL, Nagoya University

BASE_NAME=WRF_proc_AD
IST=1   # noise start
IEN=100 # noise end

declare -a bgnoiseRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")

for bgnoiseRate in "${bgnoiseRates[@]}"; do
    echo "=== Gathering outputs for bgnoiseRate = ${bgnoiseRate} ==="

    output_dir="../outputdir/woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}"
    mkdir -p "$output_dir"

    for ((i = IST; i <= IEN; i++)); do
        src="./woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/${BASE_NAME}${i}/wrfout_d01_2018-07-05_12:00:00"
        dest="${output_dir}/wrfout_d01_2018-07-05_120000_woinput_AD${i}"

        if [ -f "$src" ]; then
            cp "$src" "$dest"
            echo "case no ${i} gathered"
        else
            echo "wrfout not found for case ${i} at bgnoiseRate=${bgnoiseRate}"
        fi
    done
done
