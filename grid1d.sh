#!/bin/bash

#SBATCH --job-name=HK_grid1d
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --array=0-2
#SBATCH --output=outputs/logs/grid1d_%A_%a.out
#SBATCH --mem-per-cpu=4G
#SBATCH --time=8:00:00

# Submit from this directory: sbatch grid1d.sh

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/grid1d/

PAR=(dcp sin23 dm32)

P=${PAR[$SLURM_ARRAY_TASK_ID]}

gundamFitter -c config.yaml -of ./inputs/fitter/useGrid1D_$P.yaml -O /fitterEngineConfig/generateSamplePlots=false -o outputs/grid1d/grid1d_$P.root -t 8
