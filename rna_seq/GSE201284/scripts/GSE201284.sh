#!/bin/bash
#SBATCH --job-name=GSE201284
#SBATCH --partition=amd-hdr100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=33G
#SBATCH --time=6-06:00:00
#SBATCH --output=/home/apatil3/cancer/datasets/GSE201284/scripts/GSE201284_%j.out
#SBATCH --error=/home/apatil3/cancer/datasets/GSE201284/scripts/GSE201284_%j.err

set -euo pipefail

PROJECT_ROOT="/home/apatil3/cancer/datasets/GSE201284"
PYTHON="/home/apatil3/miniconda3/envs/genecorr/bin/python"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PYTHONUNBUFFERED=1

cd "${PROJECT_ROOT}"

"${PYTHON}" \
  "${PROJECT_ROOT}/scripts/GSE201284_all_feature_correlations.py" \
  --tpm "${PROJECT_ROOT}/raw_files/GSE201284_processed_TPM.txt" \
  --metadata "${PROJECT_ROOT}/raw_files/sample_metadata.tsv" \
  --output-root "${PROJECT_ROOT}/GSE201284_all_results" \
  --chunk-size 128
