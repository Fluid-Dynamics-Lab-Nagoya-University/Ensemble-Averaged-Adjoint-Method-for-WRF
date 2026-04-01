#!/bin/bash
#$ -P FS01OCT22
#$ -jc single
#$ -N WRF_maware
#$ -cwd
#$ -V

. /etc/profile.d/modules.sh
module load intel/2023.0.0 gcc/8.2.0
#module load intel/2022.1.2
#echo "WRFDA started..." 
cat hpcrun_wrfarw_single.sh > log.wrfarw.$JOB_ID
echo "//////////////////// namelist.input //////////////////// " >> log.wrfarw.$JOB_ID
cat namelist.input >> log.wrfarw.$JOB_ID

echo "//////////////////// log.wrfarw //////////////////// " >> log.wrfarw.$JOB_ID
./wrf.exe >> log.wrfarw.$JOB_ID

echo "//////////////////// rsl.out //////////////////// " >> log.wrfarw.$JOB_ID
cat rsl.out.0000 >> log.wrfarw.$JOBID
mv rsl.error.0000 rsl.error.0000.$JOB_ID
mv rsl.out.0000 rsl.out.0000.$JOB_ID
echo "calculation terminated" 


