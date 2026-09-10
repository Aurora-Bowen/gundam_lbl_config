#!/bin/bash

#SBATCH --job-name=HK_grid
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --array=0-50%5
#SBATCH --output=outputs/logs/grid_%A_%a.out
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00

# Submit from this directory: sbatch grid.sh

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/grid/

gundamFitter -c config.yaml -of ./inputs/fitter/useGridSearch.yaml -O /fitterEngineConfig/generateSamplePlots=false -O /fitterEngineConfig/minimizerConfig/firstNode=$((SLURM_ARRAY_TASK_ID * 1968)) -O /fitterEngineConfig/minimizerConfig/nNodes=1968 -o outputs/grid/grid_$SLURM_ARRAY_TASK_ID.root -t 8
