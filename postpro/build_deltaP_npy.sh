#!/bin/bash
set -euo pipefail

module load netcdf-fortran
module load netcdf-c

cd "$(dirname "$0")"

FC=mpifrt
# Some Fujitsu compilers require explicit free-form declaration; -Free ensures compatibility
FFLAGS="-O2 -w -Free -Kopenmp"
LIBS="-lnetcdff -lnetcdf"

echo "[build] compiling deltaP_mat_from_wrfout ..."
$FC $FFLAGS deltaP_mat_from_wrfout.f90 $LIBS -o deltaP_mat_from_wrfout
echo "[build] done."
ls -l deltaP_mat_from_wrfout

