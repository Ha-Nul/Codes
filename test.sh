#!/bin/sh
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --output sample.o0
#SBATCH --error sample.e0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skyshn95@dgist.ac.kr

module run gcc

./code

exit 0

