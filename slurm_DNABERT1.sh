#!/bin/bash

# Job name:
#SBATCH --job-name=dnabert2_embed

# Request nodes:
#SBATCH --nodelist=node08

# Processors per task:
#SBATCH --ntasks=1

# CPUs per task:
#SBATCH --cpus-per-task=2

# Memory for the job
#SBATCH --mem=10G

# Wall clock limit:
#SBATCH --time=200:00:00

# Partition:
#SBATCH --partition=cpu

# Array: 216 jobs (0-215), max 5 running at once
#SBATCH --array=0-215%5

# Output files:
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/dnabert2-%A_%a.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/dnabert2-%A_%a.err

set -euo pipefail

# -------------------------------------------------------------------------------------- #
# Environment
# -------------------------------------------------------------------------------------- #

source /home/ctools/opt/anaconda-2025-12-2/bin/activate
conda activate /net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM

# -------------------------------------------------------------------------------------- #
# Paths
# -------------------------------------------------------------------------------------- #

PROJECT_DIR="/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction"
SCRIPT_PATH="${PROJECT_DIR}/DNABERT1_embeddings.py"
DATASET_LIST="${PROJECT_DIR}/all_studies.txt"

cd "${PROJECT_DIR}"

DATASET=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${DATASET_LIST}")

if [[ -z "${DATASET}" ]]; then
    echo "No dataset found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# -------------------------------------------------------------------------------------- #
# Job
# -------------------------------------------------------------------------------------- #

echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "DATASET=${DATASET}"
echo "START=$(date)"

/usr/bin/time -v python3 "${SCRIPT_PATH}" "${DATASET}"

echo "END=$(date)"
