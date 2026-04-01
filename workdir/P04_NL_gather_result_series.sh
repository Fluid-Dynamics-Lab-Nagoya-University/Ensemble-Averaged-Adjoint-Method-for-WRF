#!/bin/bash
# Collects WRF nonlinear output files from each ensemble member's working directory
# and gathers them into a single output directory per parameter combination.
#
# by Shan Jiang, FDL, Nagoya University

BASE_NAME=WRF_proc_NL
IST=1   # noise case start

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
        output_dir="../outputdir/absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}"
        mkdir -p "$output_dir"

        for ((i = IST; i <= IEN; i++)); do
            src="./absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}/${BASE_NAME}${i}/wrfout_d01_2018-07-05_12:00:00"
            dest="${output_dir}/wrfout_d01_2018-07-05_120000_woinput_NL${i}"

            cp "$src" "$dest"
            echo "case no ${i} gathered for ng=${noiseValidRate}, ig=${inputRate}${ig_suffix}, sg=${noiseSensRate}"
        done
    done
done
