#!/bin/bash
# Distributes perturbed WRF initial conditions and adjoint forcing files
# into each ensemble member's working directory on HPC.
#
# by Shan Jiang, FDL, Nagoya University

BASE_NAME=WRF_proc_AD
IST=1   # noise case start
IEN=100 # noise case end

declare -a bgnoiseRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")

for bgnoiseRate in "${bgnoiseRates[@]}"; do
    echo "=== Processing bgnoiseRate = ${bgnoiseRate} ==="

    for ((i = IST; i <= IEN; i++)); do
        src_input="../inputdir/woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/wrfinput_d01_woinput_${i}"
        dst_input="./woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/${BASE_NAME}${i}/wrfinput_d01"

        src_rainnc="../prepare/woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/A_RAINNC_absG_i0M.dat"
        dst_rainnc="./woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/${BASE_NAME}${i}/iodir/A_RAINNC.dat"

        rsync "$src_input" "$dst_input"
        rsync "$src_rainnc" "$dst_rainnc"

        echo "copied case no $i for bgnoiseRate = ${bgnoiseRate}"
    done
done
