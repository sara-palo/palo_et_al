#! /bin/bash
#SBATCH -n 4 
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH -A lu2026-2-34
#SBATCH -J padding
#SBATCH -o padding.%j.out
#SBATCH -e padding.%j.err
module purge

module load Anaconda3/2024.06-1
conda activate scenicplus 

module load  GCC/12.3.0 BEDTools/2.31.0

OUT_DIR="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/merged"
REGION_BED="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/merged/mergedpeaks.bed"
GENOME_FASTA="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/fasta/hg38.fa"
CHROMSIZES="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/fasta/hg38.chrom.sizes"
DATABASE_PREFIX="merged_1kb"
SCRIPT_DIR=""${PWD}"/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${OUT_DIR}/merged_1kb_padding.fa \
        1000 \
        yes
