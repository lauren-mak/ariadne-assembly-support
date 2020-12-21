#!/bin/bash
#SBATCH --partition=panda					# IH lab partition
#SBATCH --nodes=1
#SBATCH --mem=100gb							# Job memory request
#SBATCH --ntasks=8

# $1 = sample
# $2 = type
# $3 = chromosome
HOME='/athena/ihlab/scratch/lam4003'
PREFIX="${1}-${3}_${2}"
READIR="${HOME}/${1}_${2}_reads"

spack load samtools/qr4zqdd
# samtools view -hb /athena/ihlab/scratch/dmm2017/chm1-13/CHM1_CHM13/CHM1_180GB_CrG_GRCh38_phased_possorted.bam "chr${3}" > ${READIR}/${PREFIX}.bam

# Converts linked-read BAM files to a mkfastq-compatible directory of FastQs. There are R1/R2 paired-end files, but the barcodes are in the I1 files and thus unattached. 
# Made with PN:longranger.lariat but already has BAM to FastQ tags, for some reason?
# /home/lam4003/bin/bamtofastq-1.2.0 --nthreads=8 ${READIR}/${PREFIX}.bam "${READIR}/chr${3}_raw"

# Append the barcode (under the BAM TAGS column, now in the I1 files) to the read name.
# Separate reads are in chm1_10x_reads/chr21_raw/40589_MissingLibrary_1_unknown_fc but longranger can find it. Need to specify full path for 'id'.
/home/dmm2017/longranger-2.1.5/longranger basic --localcores=8 --fastqs="${READIR}/chr${3}_raw" --id="chr${3}_dataset"

mv "chr${3}_dataset" "${READIR}"

# Unzip and deinterleave output FastQ file from Long Ranger. 
gunzip "${READIR}/chr${3}_dataset/outs/barcoded.fastq.gz"
${HOME}/scripts/deinterleave_fastq.sh < "${READIR}/chr${3}_dataset/outs/barcoded.fastq" ${READIR}/${PREFIX}.R1.fastq ${READIR}/${PREFIX}.R2.fastq