#!/bin/bash
#SBATCH --job-name=wrapup
#SBATCH --nodelist=node08
#SBATCH --output=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/wrapup-%j.out
#SBATCH --error=/home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction/logs/wrapup-%j.err
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G

set -euo pipefail

cd /home/projects2/kvs_students/2026/jl_qtl_prediction/repo/qtl_prediction

/net/well/pool/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python 00run1.py
