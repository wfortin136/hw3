#!/bin/bash
#SBATCH --job-name=openmp_analysis

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=hw4_analysis.out

export export OMP_NUM_THREADS=1
./serial 9000000 1000
./atomic 9000000 1000
./tpriv 9000000 1000

export export OMP_NUM_THREADS=2
./serial 9000000 1000
./atomic 9000000 1000
./tpriv 9000000 100

export export OMP_NUM_THREADS=4
./serial 9000000 1000
./atomic 9000000 1000
./tpriv 9000000 100

export export OMP_NUM_THREADS=8
./serial 9000000 1000
./atomic 9000000 1000
./tpriv 9000000 100

export export OMP_NUM_THREADS=16
./serial 9000000 1000
./atomic 9000000 1000
./tpriv 9000000 100

