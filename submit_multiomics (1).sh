#!/bin/bash
#SBATCH --job-name=cfMultiOmics
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=120:00:00
#SBATCH --partition=general

set -euo pipefail

# ── User-configurable paths ────────────────────────────────────────────────────
# Set PROJECT_DIR to the absolute path of your cloned repository on the HPC.
# Set CONDA_SH to the conda.sh initialisation script for your conda installation.
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
# ──────────────────────────────────────────────────────────────────────────────

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

# Load Nextflow environment
source "${CONDA_SH}"
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
