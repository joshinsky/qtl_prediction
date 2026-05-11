#!/bin/bash

# Job name:
#SBATCH --job-name=qtl_classifier_final_eval

# Partition & Resources
#SBATCH --partition=gpu
#SBATCH --gres=shard:1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Output files:
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/classifier/eval.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/classifier/eval.err

set -euo pipefail

# -------------------------------------------------------------------------------------- #
# Environment & Paths
# -------------------------------------------------------------------------------------- #

source /home/ctools/opt/anaconda-2025-12-2/bin/activate
conda activate /net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM

PROJECT_DIR="/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction"
SCRIPT_PATH="${PROJECT_DIR}/src/classifier/classifier.py"

cd "${PROJECT_DIR}"

# Set the target folder path
OUT_BASE="xgboost_final_eval"
RESULTS_DIR="results/output/classifier/${OUT_BASE}"
FIGURES_DIR="results/figures/${OUT_BASE}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${FIGURES_DIR}"

# -------------------------------------------------------------------------------------- #
# Job Execution
# -------------------------------------------------------------------------------------- #

/usr/bin/time -v python3 "${SCRIPT_PATH}" \
    --classifier "xgboost" \
    --pca "auto" \
    --window_size "1000" \
    --embedding_type "alt" \
    --target_label "standard" \
    --class_weighting "weighted" \
    --eval_set "test" \
    --outfile "${OUT_BASE}/${OUT_BASE}"

