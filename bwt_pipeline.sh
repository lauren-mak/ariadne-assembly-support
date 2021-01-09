#!/bin/bash
#SBATCH --partition=panda

luigid --background --logdir /athena/masonlab/scratch/users/lam4003/ariadne_data/luigi_logs
PYTHONPATH='/home/lam4003/bin/scripts/ariadne_assembly_support' luigi --module bwt_pipeline bwt_de_Novo_Assembly --workers ${1}
