#!/bin/bash
#SBATCH --job-name=gfdl_rf24_3x3
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-gfdl-ito2024-rf-3x3-%j.out
#SBATCH --error=slurm-gfdl-ito2024-rf-3x3-%j.err

set -euo pipefail
cd /scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures_p2
export MPLCONFIGDIR=/tmp/mpl-gfdl-ito2024-rf-${SLURM_JOB_ID}
mkdir -p "$MPLCONFIGDIR"
/home/GEOCLIM/software/anaconda3/2024.10/envs/geoclim/bin/python make_gfdl_ito2024_rf_3x3.py
