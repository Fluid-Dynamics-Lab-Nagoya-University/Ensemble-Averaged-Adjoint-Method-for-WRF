#!/bin/bash
# Generates and submits WRF nonlinear job scripts to the Fujitsu HPC scheduler (pjsub).
# Each job covers a batch of ensemble members for one noise rate and input rate combination.
#
# by Shan Jiang, FDL, Nagoya University

cases_per_job=20

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
        total_cases=1
    else
        total_cases=100
    fi

    for inputRate in "${inputRates[@]}"; do
        for ((start=1; start<=total_cases; start+=cases_per_job)); do
            end=$((start + cases_per_job - 1))
            ((end > total_cases)) && end=$total_cases

            job_name="job_NL_absG_i0M_NOV0.02_ng${noiseValidRate}_ig${inputRate}${ig_suffix}_sg${noiseSensRate}_${start}-${end}.sh"

            cat << EOF > $job_name
#!/bin/bash
#PJM -L rscgrp=fx-small
#PJM -L node=3
#PJM --mpi proc=96
#PJM -L elapse=48:00:00
#PJM -j
#PJM -N NL${start}-${end}

module load wrf/3.9.1.1

noiseSensRate=${noiseSensRate}
noiseValidRate=${noiseValidRate}
inputRate=${inputRate}
ig_suffix=${ig_suffix}

BASE_NAME=WRF_proc_NL
IST=${start}
IEN=${end}

for ((i = \$IST; i <= \$IEN; i++)); do
  cd absG_i0M_NOV0.02_ng\${noiseValidRate}_ig\${inputRate}\${ig_suffix}_sg\${noiseSensRate}/\${BASE_NAME}\${i}
  mpiexec -n 96 ./wrf.exe
  cd - > /dev/null
  echo "finished calc for case no \${i}"
done
EOF

            # submit job
            pjsub $job_name
        done
    done
done
