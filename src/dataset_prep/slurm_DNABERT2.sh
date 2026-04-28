#!/bin/bash

# Job name:
#SBATCH --job-name=dnabert2_embed

# Partition:
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1 

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

/usr/bin/time -v python3 "${SCRIPT_PATH}" "${DATA_PATH}"
