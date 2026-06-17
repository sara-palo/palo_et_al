#! /bin/bash
#SBATCH -n 48 
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -A lu2026-2-34
#SBATCH -J cistarget
#SBATCH -o cistarget.%j.out
#SBATCH -e cistarget.%j.err
module purge

module load Anaconda3/2024.06-1
conda activate scenicplus 

WORK_DIR=""${PWD}""
OUT_DIR="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/cistarget_outs"
CBDIR="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/aertslab_motif_collection/v10nr_clust_public/singletons"
FASTA_FILE="/home/spalo/lu2025-12-50/spalo/scenicplus_storage/merged/merged_1kb_padding.fa"
MOTIF_LIST="${WORK_DIR}/motifs.txt"
SCRIPT_DIR="${WORK_DIR}/create_cisTarget_databases"
DATABASE_PREFIX="cisTarget_database_merged"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20

