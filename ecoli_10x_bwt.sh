#!/bin/bash
#SBATCH --partition=panda
#SBATCH --job-name=ecoli_10x_bwt
#SBATCH --output=ecoli_10x_bwt.out
#SBATCH --mem=20gb
#SBATCH --ntasks=20

luigid --background
PYTHONPATH='.' luigi --module /home/lam4003/bin/scripts/ariadne_assembly_support/bwt_pipeline bwt_de_Novo_Assembly --workers 20
