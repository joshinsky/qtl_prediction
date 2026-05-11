#!/bin/bash

# Job name:
#SBATCH --job-name=qtl_classifier

# Partition & Resources
#SBATCH --partition=gpu
#SBATCH --gres=shard:1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Job array 
#SBATCH --array=0-15%8

# Output files:
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/classifier/class-%A_%a.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/classifier/class-%A_%a.err

set -euo pipefail

# -------------------------------------------------------------------------------------- #
# Environment & Paths
# -------------------------------------------------------------------------------------- #

source /home/ctools/opt/anaconda-2025-12-2/bin/activate
conda activate /net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/repo/dnaLM/conda_env/dnaLM

PROJECT_DIR="/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction"
SCRIPT_PATH="${PROJECT_DIR}/src/classifier/classifier.py"

cd "${PROJECT_DIR}"

# -------------------------------------------------------------------------------------- #
# Parameter Definitions
# -------------------------------------------------------------------------------------- #

# Define the arrays for each parameter
CLASSIFIERS=("xgboost")
PCAS=("skip" "auto")
WINDOWS=("20" "100")
EMBEDDINGS=("alt")
TARGETS=("standard" "single")
WEIGHTINGS=("weighted" "none")

# -------------------------------------------------------------------------------------- #
# Map Task ID to Parameters using modulo arithmetic
# -------------------------------------------------------------------------------------- #

# Create an index mapping
IDX=$SLURM_ARRAY_TASK_ID

# Extract parameter based on modulus and division

w_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
W="${WEIGHTINGS[$w_idx]}"

t_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
T="${TARGETS[$t_idx]}"

emb_idx=$(( IDX % 1 )); IDX=$(( IDX / 1 ))
EMB="${EMBEDDINGS[$emb_idx]}"

win_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
WIN="${WINDOWS[$win_idx]}"

pca_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
PCA="${PCAS[$pca_idx]}"

# Define output file name/folder name
OUT_BASE="xgboost_wt-${W}_tgt-${T}_pca-${PCA}_win-${WIN}_emb-${EMB}"

# Set the target folder path
RESULTS_DIR="results/output/classifier/${OUT_BASE}"
FIGURES_DIR="results/figures/${OUT_BASE}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${FIGURES_DIR}"

echo "Running task $SLURM_ARRAY_TASK_ID with parameters:"
echo "Classifier: xgboost | PCA: $PCA | Window: $WIN | Embedding: $EMB | Target: $T | Weighting: $W"
echo "Output file: $OUT_BASE"

# -------------------------------------------------------------------------------------- #
# Job Execution
# -------------------------------------------------------------------------------------- #

/usr/bin/time -v python3 "${SCRIPT_PATH}" \
    --classifier "xgboost" \
    --pca "$PCA" \
    --window_size "$WIN" \
    --embedding_type "$EMB" \
    --target_label "$T" \
    --class_weighting "$W" \
    --eval_set "val" \
    --outfile "${OUT_BASE}/${OUT_BASE}"

