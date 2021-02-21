#!/bin/bash
#SBATCH --partition=panda

# $1 = assembly mode (ariadne, bowtie)
# $2 = Luigi config path
# $3 = number of threads (only works on curie)

spack load gcc@6.3.0
export LUIGI_CONFIG_PATH=${2}
luigid --background --logdir /athena/masonlab/scratch/users/lam4003/ariadne_data/luigi_logs
PYTHONPATH='/home/lam4003/bin/scripts/ariadne_assembly_support' luigi --module ${1}_pipeline de_Novo_Assembly --workers ${3}
