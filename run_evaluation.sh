#!/bin/bash
#SBATCH --partition=panda

# Set basics
# $1 = dataset prefix 
# $2 = master directory
# $3 = reference FastA or directory
ANALYSES_DIR="${2}/${1}_analyses"
OUTPUT_DIR="${2}/${1}_output"

# metaQUAST scaffold comparison
python /home/lam4003/bin/quast-master/metaquast.py \
    -l 'Illumina,Species,Fragments,5000_Old,5000,10000,15000,20000' \
    -o ${ANALYSES_DIR}/metaQUAST -r ${3} \
    ${ANALYSES_DIR}/illumina.scaffolds.fasta \
    ${ANALYSES_DIR}/spc.scaffolds.fasta \
    ${ANALYSES_DIR}/frg.scaffolds.fasta \
    ${ANALYSES_DIR}/5000_old.scaffolds.fasta \
    ${ANALYSES_DIR}/5000.scaffolds.fasta \
    ${ANALYSES_DIR}/10000.scaffolds.fasta \
    ${ANALYSES_DIR}/15000.scaffolds.fasta \
    ${ANALYSES_DIR}/20000.scaffolds.fasta
    # -l No_Deconv
    # ${ANALYSES_DIR}/no_deconv.scaffolds.fasta \

# Read cloud summary statistics
python3 /home/lam4003/bin/scripts/ariadne_assembly_support/cli.py evaluate_clouds \
    -d 'Species,Fragments,5000,10000,15000,20000' \
    -p ${1}_spc_bsort,${1}_frg_bsort,5000,10000,15000,20000 \
    ${ANALYSES_DIR}/
    # -d No_Deconv
    # -p no_deconv

cp ${ANALYSES_DIR}/*Purity* ${OUTPUT_DIR}/
cp ${ANALYSES_DIR}/*Entropy* ${OUTPUT_DIR}/
cp ${ANALYSES_DIR}/Avg_Cloud_Size.png ${OUTPUT_DIR}/
cp ${ANALYSES_DIR}/Avg_Summary_Stats.tbl ${OUTPUT_DIR}/
cp ${ANALYSES_DIR}/metaQUAST/combined_reference/report.tsv ${OUTPUT_DIR}/
cp ${ANALYSES_DIR}/metaQUAST/*html ${OUTPUT_DIR}/
cp -r ${ANALYSES_DIR}/metaQUAST/icarus_viewers ${OUTPUT_DIR}/
