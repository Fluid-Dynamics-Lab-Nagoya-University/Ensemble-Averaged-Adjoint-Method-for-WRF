#!/bin/bash
# Re-links wrf.exe into each ensemble member's working directory for NL runs.
# Edit your_wrf_directory to point to the parent folder of your local WRF installation.
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
        for ((i = IST; i <= IEN; i++)); do
            target_dir="./absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}/${BASE_NAME}${i}"

            if [ -d "$target_dir" ]; then
                cd "$target_dir" || { echo "cd failed: $target_dir"; exit 1; }
                rm -f wrf.exe
                ln -s your_wrf_directory/WRF/WRFPLUS-3.9.1.1_mdfv2/main/wrf.exe
                cd - > /dev/null
            else
                echo "Directory not found: $target_dir"
            fi
        done
    done
done
