#!/bin/bash
#SBATCH -J multimethod_full_o2
#SBATCH --output=multimethod_full_o2_%j.out
#SBATCH --error=multimethod_full_o2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --account=lrgroup
set -euo pipefail
source /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/etc/profile.d/conda.sh
conda activate /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/envs/xmip
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/tmp/matplotlib-db9274
cd /scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures_p2
python build_multimethod_mean_full_column.py
