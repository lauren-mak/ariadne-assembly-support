#!/bin/bash
#SBATCH --partition=panda                   # IH lab partition
#SBATCH --ntasks=30
#SBATCH --mem=150gb                         # Job memory request

# $1 = sample
# $2 = type
# $3 = search distance
HOME='/athena/masonlab/scratch/users/lam4003/ariadne_data'
PREFIX="${1}_${2}"
OUTDIR="${HOME}/${PREFIX}_${3}"
IHDIR='/athena/ihlab/scratch/lam4003'
READIR="${IHDIR}/mock_original_reads/${PREFIX}"

spack load gcc@6.3.0

echo "Continuing cloudSPAdes + Ariadne at k = 55 for the ${1} ${2} dataset with a search distance of ${3}..."
/home/lam4003/bin/spades/assembler/spades.py --restart-from k55 --search-distance ${3} -o ${OUTDIR}

mv ${OUTDIR}/K55/"${3}".* ${OUTDIR}/

echo "Completing the ${1} ${2} with search distance ${3} dataset..."
python ${IHDIR}/scripts/ariadne_support.py complete_reads ${READIR}.fixed.R1.fastq ${OUTDIR}/${3}.R1.fastq ${OUTDIR}/${3}_full.R1.fastq
python ${IHDIR}/scripts/ariadne_support.py complete_reads ${READIR}.fixed.R2.fastq ${OUTDIR}/${3}.R2.fastq ${OUTDIR}/${3}_full.R2.fastq

sed -i "s|K55/${3}|${3}_full|g" ${OUTDIR}/input_dataset.yaml

echo "Continuing cloudSPAdes + Ariadne for the ${1} ${2} dataset with a search distance of ${3}..."
/home/lam4003/bin/spades/assembler/spades.py --restart-from k55 --search-distance 0 -o ${OUTDIR}
