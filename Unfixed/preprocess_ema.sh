#!/bin/bash
#SBATCH --partition=panda                   # IH lab partition
"#SBATCH --ntasks=10
#SBATCH --mem=200gb                         # Job memory request"

# $1 = sample
# $2 = type
HOME='/athena/ihlab/scratch/lam4003'
PREFIX="${1}_${2}"
READIR="${HOME}/mock_original_reads/${PREFIX}"
WHITELIST="${HOME}/chm1_10x_reads/10x_bc_whitelist.txt"
OUTDIR="/athena/masonlab/scratch/users/lam4003/ariadne_data/${PREFIX}_ema"
OUT_PREFIX="${OUTDIR}/${PREFIX}"

# EMA_CMD="/home/lam4003/bin/ema/ema"
# mkdir $OUTDIR

# echo "Interleaving reads from the ${1} ${2} dataset..."
# paste ${READIR}.R1.fastq ${READIR}.R2.fastq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' > ${OUT_PREFIX}.interleave.fastq

# echo "Generating EMA count files for the ${1} ${2} dataset..."
# cat ${OUT_PREFIX}.interleave.fastq | ${EMA_CMD} count -w ${WHITELIST} -o ${OUT_PREFIX}

# echo "Generating EMA read-bins for the ${1} ${2} dataset..."
# cat ${OUT_PREFIX}.interleave.fastq | ${EMA_CMD} preproc -w ${WHITELIST} -o ${OUTDIR} ${OUT_PREFIX}.ema-ncnt

echo "Running EMA for the ${1} ${2} dataset with 500 bins..."
LEN_TRUE=3
for i in $(seq 0 499)
do
    BIN_NUM="${i}"
    LEN_I=${#BIN_NUM}
    LEN_SHORT=$(expr $LEN_TRUE - $LEN_I)
    if [ LEN_SHORT > 0 ]
    then
        for j in $(seq 1 $LEN_SHORT)
        do
            BIN_NUM="0${BIN_NUM}"
        done
    fi
    sbatch -J ema_${BIN_NUM} -o ${OUTDIR}/${BIN_NUM}.out --mem=1gb ${HOME}/scripts/deconvolve_ema.sh ${1} ${OUTDIR}/ema-bin-${BIN_NUM}
done
