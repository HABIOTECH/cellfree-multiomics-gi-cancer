#!/bin/bash
#SBATCH --job-name=cfMultiOmics
#SBATCH --output=/sc/kzd/home/abdullah/bioinformatics/CRC-NC-multiomics/logs/pipeline_%j.out
#SBATCH --error=/sc/kzd/home/abdullah/bioinformatics/CRC-NC-multiomics/logs/pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=120:00:00
#SBATCH --partition=general

set -euo pipefail

# Project Directory
PROJECT_DIR="/sc/kzd/home/abdullah/bioinformatics/CRC-NC-multiomics"
mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

# Load Nextflow environment
source /sc/kzd/home/abdullah/miniconda/etc/profile.d/conda.sh
conda activate nxf_env

# Run Nextflow
nextflow run main.nf \
    -c nextflow.config \
    -profile cluster \
    -resume \
    --dna_csv SraRunTable_cfDNA.csv \
    --medip_csv SraRunTable_cfMeDIP.csv \
    --total_rna_csv SraRunTable_SRR_Acc_List_total_cfRNA.csv \
    --small_rna_csv SraRunTable_Small_cfRNA.csv \
    -work-dir "${PROJECT_DIR}/work"
