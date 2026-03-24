#!/bin/bash

# Job name:
#SBATCH --job-name=qtl_pipeline_run1_fix

# Request nodes:
#SBATCH --nodes=1

# Processors per task:
#SBATCH --ntasks=1

# CPUs per task:
#SBATCH --cpus-per-task=1

# Memory for the job
#SBATCH --mem=20G

# Wall clock limit
#SBATCH --time=200:00:00

# Partition:
#SBATCH --partition=cpu

#SBATCH --nodelist=node08

# Output files:
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/run1_fix-%j.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/run1_fix-%j.err

# -------------------------------------------------------------------------------------- #
# Environment
# -------------------------------------------------------------------------------------- #

source /home/ctools/opt/anaconda-2025-12-2/bin/activate
conda activate /net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env

# -------------------------------------------------------------------------------------- #
# Job Execution
# -------------------------------------------------------------------------------------- #

# Navigate to the directory containing the run1.py script
cd /home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction

echo "Re-running 00run1.py after filename fix..."

# Run the corresponding Python script based on the array ID
/usr/bin/time -v /net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3 00run1.py

echo "Finished array task ${SLURM_ARRAY_TASK_ID}."
