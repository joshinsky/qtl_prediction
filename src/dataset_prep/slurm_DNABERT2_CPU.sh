#!/bin/bash

# Job name:
#SBATCH --job-name=dnabert2_embed_cpu

# Partition:
#SBATCH --partition=cpu

# Array configuration for 3 tasks:
#SBATCH --array=0-2

# Processors per task:
#SBATCH --ntasks=1

# CPUs per task (increased to 8 to speed up CPU execution):
#SBATCH --cpus-per-task=10

# Memory for the job
#SBATCH --mem=32G

# Wall clock limit:
#SBATCH --time=200:00:00

# Output files:
#SBATCH --output=logs/dnabert2-cpu-%A_%a.out
#SBATCH --error=logs/dnabert2-cpu-%A_%a.err

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
SCRIPT_PATH="${PROJECT_DIR}/src/dataset_prep/DNABERT2_embeddings_CPU.py"
DATA_PATH="${PROJECT_DIR}/results/output/dataset_prep/final_full_dataset.tsv.gz"

cd "${PROJECT_DIR}"

# -------------------------------------------------------------------------------------- #
# Job
# -------------------------------------------------------------------------------------- #

WINDOWS=(20 1000 100)
CURRENT_WINDOW="${WINDOWS[$SLURM_ARRAY_TASK_ID]}"

# Tell PyTorch to fully utilize the 8 requested CPU cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

/usr/bin/time -v python3 "${SCRIPT_PATH}" "${DATA_PATH}" "${CURRENT_WINDOW}"
