#!/bin/bash

# Job name:
#SBATCH --job-name=train_test_split

# Partition:
#SBATCH --partition=cpu

# Processors per task:
#SBATCH --ntasks=1

# CPUs per task:
#SBATCH --cpus-per-task=4

# Memory for the job
#SBATCH --mem=32G

# Wall clock limit:
#SBATCH --time=200:00:00

# Output files:
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/split-%A_%a.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/split-%A_%a.err

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
SCRIPT_PATH="${PROJECT_DIR}/src/dataset_prep/holdout_testset.py"

cd "${PROJECT_DIR}"

# -------------------------------------------------------------------------------------- #
# Job
# -------------------------------------------------------------------------------------- #

# Run the python script
/usr/bin/time -v python3 "${SCRIPT_PATH}"
