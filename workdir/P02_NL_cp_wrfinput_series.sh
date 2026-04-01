#!/bin/bash
# Distributes perturbed WRF initial conditions into each ensemble member's
# working directory on HPC for NL runs.
#
# by Shan Jiang, FDL, Nagoya University

BASE_NAME=WRF_proc_NL
IST=1   # noise index start

# input_mode controls the suffix appended after ig${inputRate} in directory names:
#   0: no suffix        (both positive and negative input)
#   1: suffix "P"       (positive input only)
#   2: suffix "N"       (negative input only)
#   3: suffix "Nactr"   (negative input at the 8 selected actuator points only)
input_mode=1

case "$input_mode" in
    0) ig_suffix="" ;;
    1) ig_suffix="P" ;;
    2) ig_suffix="N" ;;
    3) ig_suffix="Nactr" ;;
    *) echo "Error: invalid input_mode=${input_mode}"; exit 1 ;;
esac

declare -a noiseValidRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")
declare -a inputRates=("0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9")

for rate in "${noiseValidRates[@]}"; do
    noiseSensRate="$rate"
    noiseValidRate="$rate"

    # Determine ensemble size based on noise rate
    if [[ "$noiseValidRate" == "0" ]]; then
        IEN=1
    else
        IEN=100
    fi

    for inputRate in "${inputRates[@]}"; do
        indir="../inputdir/absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}"
        workdir="./absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}"

        for ((i = IST; i <= IEN; i++)); do
            src="${indir}/wrfinput_d01_woinput_${i}"
            dest="${workdir}/${BASE_NAME}${i}/wrfinput_d01"

            rsync "$src" "$dest"
            echo "copied case no $i for ng=${noiseValidRate}, ig=${inputRate}${ig_suffix}, sg=${noiseSensRate}"
        done
    done
done
