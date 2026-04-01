#!/bin/bash
# Re-links wrf.exe into each ensemble member's working directory.
# Edit your_wrf_directory to point to the parent folder of your local WRF installation.
#
# by Shan Jiang, FDL, Nagoya University

IST=1   # noise case start
IEN=100 # noise case end

declare -a bgnoiseRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")

for bgnoiseRate in "${bgnoiseRates[@]}"; do
    echo "=== Linking wrf.exe for bgnoiseRate = ${bgnoiseRate} ==="

    for ((i = IST; i <= IEN; i++)); do
        target_dir="woinput_absG_i0M_NOV0.02_sg_${bgnoiseRate}/WRF_proc_AD${i}"

        if [ -d "$target_dir" ]; then
            echo "Processing: $target_dir"
            cd "$target_dir" || { echo "cd failed: $target_dir"; exit 1; }

            rm -f wrf.exe
            ln -s your_wrf_directory/WRFPLUS-3.9.1.1_mdfv2/main/wrf.exe

            cd - > /dev/null
        else
            echo "Directory not found: $target_dir"
        fi
    done
done
