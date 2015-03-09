#!/bin/bash
#SBATCH --job-name=openmp_analysis

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=hw4_analysis.out

export export OMP_NUM_THREADS=1
./serial 500000 250
./atomic 500000 250
./tpriv 500000 250
./serial 250000 250
./atomic 250000 250
./tpriv 250000 250
./serial 100000 250
./atomic 100000 250
./tpriv 100000 250

export export OMP_NUM_THREADS=2
./serial 500000 250
./atomic 500000 250
./tpriv 500000 250
./serial 250000 250
./atomic 250000 250
./tpriv 250000 250
./serial 100000 250
./atomic 100000 250
./tpriv 100000 250

export export OMP_NUM_THREADS=4
./serial 500000 250
./atomic 500000 250
./tpriv 500000 250
./serial 250000 250
./atomic 250000 250
./tpriv 250000 250
./serial 100000 250
./atomic 100000 250
./tpriv 100000 250

export export OMP_NUM_THREADS=8
./serial 500000 250
./atomic 500000 250
./tpriv 500000 250
./serial 250000 250
./atomic 250000 250
./tpriv 250000 250
./serial 100000 250
./atomic 100000 250
./tpriv 100000 250

export export OMP_NUM_THREADS=16
./serial 500000 250
./atomic 500000 250
./tpriv 500000 250
./serial 250000 250
./atomic 250000 250
./tpriv 250000 250
./serial 100000 250
./atomic 100000 250
./tpriv 100000 250

