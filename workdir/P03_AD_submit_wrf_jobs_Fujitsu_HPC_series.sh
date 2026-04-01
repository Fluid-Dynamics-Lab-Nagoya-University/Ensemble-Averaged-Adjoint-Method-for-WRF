#!/bin/bash
# Generates and submits WRF adjoint job scripts to the Fujitsu HPC scheduler (pjsub).
# Each job covers a batch of ensemble members for one noise rate.
#
# by Shan Jiang, FDL, Nagoya University

total_cases=100
cases_per_job=20

BASE_NAME=WRF_proc_AD
bgnoiseRates=("0" "0.001" "0.01" "0.05" "0.1" "0.15" "0.2" "0.3" "0.4" "0.5")

for bgnoiseRate in "${bgnoiseRates[@]}"; do
  for ((start=1; start<=total_cases; start+=cases_per_job)); do
    end=$((start + cases_per_job - 1))
    ((end > total_cases)) && end=$total_cases

    # 1 script for 1 parameter combination
    script_name="job_NL_absG_i0M_NOV0.02_sg_${bgnoiseRate}_${start}-${end}.sh"

    cat << EOF > "$script_name"
#!/bin/bash
#PJM -L rscgrp=fx-small
#PJM -L node=3
#PJM --mpi proc=96
#PJM -L elapse=48:00:00
#PJM -j
#PJM -N AD${start}-${end}

module load wrf/3.9.1.1

bgnoiseRate=${bgnoiseRate}
BASE_NAME=${BASE_NAME}

IST=${start}
IEN=${end}

for ((i = \$IST; i<= \$IEN; i++)); do
  cd woinput_absG_i0M_NOV0.02_sg_\${bgnoiseRate}/\${BASE_NAME}\${i}
  mpiexec -n 96 ./wrf.exe
  cd - > /dev/null
  echo "finished calc for case no \$i"
done
EOF


    # submit
    pjsub "$script_name"
  done
done