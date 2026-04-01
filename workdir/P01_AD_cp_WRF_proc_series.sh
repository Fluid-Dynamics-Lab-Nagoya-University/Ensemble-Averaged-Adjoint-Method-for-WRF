#!/bin/bash
# Sets up base WRF working directories for HPC execution.
# Copies the WRF template directory and distributes perturbed initial conditions
# across ensemble members for each noise rate.
#
# by Shan Jiang, FDL, Nagoya University

BASE_NAME=WRF_proc_AD
IST=1   # noise member start
IEN=100 # noise member end

declare -a bgnoiseRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")

for bgnoiseRate in "${bgnoiseRates[@]}"; do
    input_dir="../inputdir/woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}"
    work_dir="./woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}"
    prepare_dir="../prepare/woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}"

    mkdir -p "$input_dir"
    mkdir -p "$work_dir"

    for ((i = IST; i <= IEN; i++)); do
        rsync -a ../"${BASE_NAME}_org/" "${work_dir}/${BASE_NAME}${i}"
        rsync "${prepare_dir}/wrfinput_d01_woinput_${i}" "${input_dir}/wrfinput_d01_woinput_${i}"
    done
done