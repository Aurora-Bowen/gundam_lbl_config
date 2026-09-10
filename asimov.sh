#!/bin/bash

#SBATCH --job-name=HK_asimov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=outputs/logs/asimov_%A.out
#SBATCH --mem-per-cpu=4G

# The plain fit, nothing profiled: all six oscillation parameters and the 70
# systematics free, minimised on the Asimov data set config.yaml defines. Its
# postfit LLH is the reference the profile and grid curves are measured from.
#
#   sbatch asimov.sh             from this directory

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/asimov/

gundamFitter -c config.yaml --out-dir outputs/asimov/ -t 8
