#!/bin/bash
#
#SBATCH --job-name=FLaMe_mpi
#SBATCH --output=res_mpi.txt
#
#SBATCH --ntasks=100
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M

module load OpenMPI/4.1.2-GCC-11.2.0
mpirun -np 99 ./runlake_paralize
