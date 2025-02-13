#!/bin/sh
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output sample.o0
#SBATCH --error sample.e0
#SBATCH --hint=nomultithread
#SBATCH --threads-per-core=1

module run gcc

./build/RSJJ_TOA test.json

exit 0 
