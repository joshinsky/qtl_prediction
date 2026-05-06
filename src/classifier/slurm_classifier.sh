#!/bin/bash

# Job name:
#SBATCH --job-name=qtl_classifier

# Partition & Resources
#SBATCH --partition=gpu
#SBATCH --gres=shard:1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=200:00:00

# Job array 
#SBATCH --array=0-31%4

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
SCRIPT_PATH="${PROJECT_DIR}/src/classifier/your_classifier_script.py"

cd "${PROJECT_DIR}"

# -------------------------------------------------------------------------------------- #
# Parameter Definitions
# -------------------------------------------------------------------------------------- #

# Define the arrays for each parameter
CLASSIFIERS=("xgboost")
PCAS=("skip" "auto")
WINDOWS=("20" "100") #"1000")
EMBEDDINGS=("alt" "delta")
TARGETS=("standard" "both")
WEIGHTINGS=("none" "weighted")

# -------------------------------------------------------------------------------------- #
# Map Task ID to Parameters using modulo arithmetic
# -------------------------------------------------------------------------------------- #

# Create an index mapping
IDX=$SLURM_ARRAY_TASK_ID

# Extract parameter based on modulus and division
# Number of choices: CLASSIFIERS(2), PCAS(2), WINDOWS(3), EMBEDDINGS(2), POSITIONS(4), TARGETS(2), WEIGHTINGS(2)

w_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
W="${WEIGHTINGS[$w_idx]}"

t_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
T="${TARGETS[$t_idx]}"

emb_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
EMB="${EMBEDDINGS[$emb_idx]}"

win_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
WIN="${WINDOWS[$win_idx]}"

pca_idx=$(( IDX % 2 )); IDX=$(( IDX / 2 ))
PCA="${PCAS[$pca_idx]}"

# Define output file name/folder name
OUT_BASE="xgboost_pca-${PCA}_win-${WIN}_emb-${EMB}_tgt-${T}_wt-${W}"

# Set the target folder path
RESULTS_DIR="${PROJECT_DIR}/results/output/classifier/${OUT_BASE}"
mkdir -p "${RESULTS_DIR}"

echo "Running task $SLURM_ARRAY_TASK_ID with parameters:"
echo "Classifier: xgboost | PCA: $PCA | Window: $WIN | Embedding: $EMB | Target: $T | Weighting: $W"
echo "Output file: $OUTFILE"

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

