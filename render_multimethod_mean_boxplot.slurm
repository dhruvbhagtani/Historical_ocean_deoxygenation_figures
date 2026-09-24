#!/bin/bash
#SBATCH -J multimethod_o2_figure
#SBATCH --output=multimethod_o2_figure_%j.out
#SBATCH --error=multimethod_o2_figure_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --account=lrgroup
set -euo pipefail
source /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/etc/profile.d/conda.sh
conda activate /scratch/gpfs/GEOCLIM/LRGROUP/db9274/Anaconda/envs/xmip
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/tmp/matplotlib-db9274
cd /scratch/gpfs/LRGROUP/db9274/Variability_trends_EC/Paper_figures_p2
python -c "import json; n=json.load(open('Oxygen_sampling_mapping_difference_trends_two_periods.ipynb')); ns={'__name__':'figure_build'}; exec(''.join(n['cells'][1]['source']),ns); exec(''.join(n['cells'][2]['source']),ns)"
