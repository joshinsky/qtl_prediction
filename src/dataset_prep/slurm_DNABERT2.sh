#!/bin/bash

# Job name:
#SBATCH --job-name=dnabert2_embed_with_ref_extended_window_range

# Partition:
#SBATCH --partition=gpu

# Request exactly 1 GPU shard (half a GPU's VRAM)
#SBATCH --gres=shard:1 

# Array configuration for 3 tasks:
#SBATCH --array=0-2

# Processors per task:
#SBATCH --ntasks=1

# CPUs per task:
#SBATCH --cpus-per-task=4

# Memory for the job
#SBATCH --mem=16G

# Wall clock limit:
#SBATCH --time=200:00:00

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
SCRIPT_PATH="${PROJECT_DIR}/src/dataset_prep/DNABERT2_embeddings.py"
DATA_PATH="${PROJECT_DIR}/results/output/dataset_prep/final_full_dataset.tsv.gz"

cd "${PROJECT_DIR}"

# -------------------------------------------------------------------------------------- #
# Job
# -------------------------------------------------------------------------------------- #

# Define the window sizes in a bash array
WINDOWS=(20 1000 100)

# Extract the current window size based on the task ID
CURRENT_WINDOW="${WINDOWS[$SLURM_ARRAY_TASK_ID]}"

/usr/bin/time -v python3 "${SCRIPT_PATH}" "${DATA_PATH}" "${CURRENT_WINDOW}"
