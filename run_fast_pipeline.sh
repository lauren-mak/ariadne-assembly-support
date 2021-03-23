#!/bin/bash
#SBATCH --partition=panda
#SBATCH --ntasks=30
#SBATCH --mem=200gb

# $1 = FastQ prefix
# $2 = search distance

PREFIX="/athena/ihlab/scratch/lam4003/microbiome_reads/${1}_bsort"

spack load gcc@6.3.0
/home/lam4003/bin/spades/assembler/spades.py --meta --only-assembler --gemcode1-1 ${PREFIX}.R1.fastq --gemcode1-2 ${PREFIX}.R2.fastq --search-distance ${2} --size-cutoff 6 -t 30 -m 200 -o /athena/masonlab/scratch/users/lam4003/ariadne_data/${1}_${2}_fast
