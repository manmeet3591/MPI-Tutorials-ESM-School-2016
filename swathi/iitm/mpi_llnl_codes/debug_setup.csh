#!/bin/tcsh -f 
#PBS -V
#PBS -l walltime=2:00:00
#PBS -l select=40:ncpus=16:mpiprocs=16
#PBS -l place=scatter:excl
#PBS -q workq

set echo on
# 05/10/10..For tutorial purposes..
cd /scratch/swathi/mpi_llnl_codes
source /usr/share/Modules/init/tcsh
module purge
module load netcdf-4.1.2
module load intel-cluster-studio-2013
module list

