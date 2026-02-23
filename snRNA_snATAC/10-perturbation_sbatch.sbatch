#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --job-name=perturb_TF
#SBATCH --output=logs/perturb_%A_%a.out
#SBATCH --error=logs/perturb_%A_%a.err
#SBATCH --time=120:00:00
#SBATCH --array=0-19%5
#SBATCH --nodelist=c3po

# inizialize conda
source /home/ciervo/miniconda3/etc/profile.d/conda.sh
conda activate scenicplus

# setwd
cd /home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

mkdir -p logs

# Load MP4 transcription factor list
# TF_LIST="/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/tf_list_MP4.txt"

TF_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TF_LIST")

python perturbation_iteration.py --tf "$TF_NAME"
