#!/bin/bash
#SBATCH --partition=panda

# $1 = assembly mode (ariadne, bowtie)
# $2 = Luigi config path
# $3 = number of threads (only works on curie)

READ_DIR="$4"
PREFIX="$5"
WORK_DIR="${PREFIX}_ema"
SORT_PREFIX="${PREFIX}_bsort"
WHITELIST="$6"
NUM_CHUNKS="$7"

# Make the working EMA directory if it doesn't exist
mkdir $WORK_DIR    

# Interleave the barcode-sorted FastQs (assuming that FastQs have been pre-sorted with some other pipeline)
python3 /home/lam4003/bin/scripts/ariadne_assembly_support/cli.py interleave_fastqs ${READ_DIR}/${SORT_PREFIX} ${WORK_DIR}

# Make EMA barcode counts
cat ${WORK_DIR}/${SORT_PREFIX}.fastq | /home/lam4003/bin/ema/ema count -w ${WHITELIST} -o ${WORK_DIR}/${SORT_PREFIX}

# Divide FastQs into EMA read-bins
cat ${WORK_DIR}/${SORT_PREFIX}.fastq | /home/lam4003/bin/ema/ema preproc -w  ${WHITELIST} -n ${NUM_CHUNKS} -o ${WORK_DIR}/bins  ${WORK_DIR}/${SORT_PREFIX}.ema-ncnt

# Rest of the EMA annotation pipeline
spack load gcc@6.3.0
export LUIGI_CONFIG_PATH=${2}
luigid --background --logdir /athena/masonlab/scratch/users/lam4003/ariadne_data/luigi_logs
PYTHONPATH='/home/lam4003/bin/scripts/ariadne_assembly_support' luigi --module ${1}_pipeline de_Novo_Assembly --workers ${3}
