#!/bin/bash
#SBATCH --partition=panda
#SBATCH --ntasks=30
#SBATCH --mem=200gb

# $1 = FastQ prefix
# $2 = original search distance
# $3 = new search distance

spack load gcc@6.3.0
/home/lam4003/bin/spades/assembler/spades.py --restart-from k55 -t 30 --search-distance ${3} -o /athena/masonlab/scratch/users/lam4003/ariadne_data/${1}_${2}/cloudSPAdes_init 