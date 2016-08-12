#!/bin/bash

#SBATCH -N 20
#SBATCH -n 160
#SBATCH -p mpidev
#SBATCH -J mexico_rap
#SBATCH -t 00:59:59


. /etc/profile.d/modules.sh
module use /cfs/devA/opt/Modules/modulefiles
module load openmpi
module load lustre
module load netcdf

mpirun -np 160 $1
