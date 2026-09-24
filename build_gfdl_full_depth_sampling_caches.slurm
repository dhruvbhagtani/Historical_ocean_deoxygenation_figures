#!/bin/bash
#SBATCH -J gfdl_full_sampling
#SBATCH --output=gfdl_full_sampling_%A_%a.out
#SBATCH --error=gfdl_full_sampling_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --account=lrgroup
set -euo pipefail
source /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/etc/profile.d/conda.sh
conda activate /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/envs/xmip
export PYTHONDONTWRITEBYTECODE=1
cd /scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures_p2
methods=(oi rf)
python build_gfdl_full_depth_sampling_caches.py --method "${methods[$SLURM_ARRAY_TASK_ID]}"
