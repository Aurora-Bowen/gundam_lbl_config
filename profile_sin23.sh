#!/bin/bash

#SBATCH --job-name=HK_prof_sin23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-49%8
#SBATCH --output=outputs/logs/profile_sin23_%A_%a.out
#SBATCH --mem-per-cpu=4G

# Profile likelihood in sin^2(theta23). Fifty independent fits, one per array
# task, each freezing sin^2(theta23) at one grid point and minimising the other
# five oscillation parameters and all 70 systematics around it.
#
# Grid: 0.300 to 0.692 in steps of 0.008, which is 0.38 sigma per step against
# the 0.021 in MakeCovariances.C. Prior 0.528, upper parameterLimit 0.700.
#
# Every task writes its own file into outputs/profile/, named from the override
# it was given, so they cannot collide. See inputs/profile/README.md.
#
#   sbatch profile_sin23.sh      from this directory

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/profile/

POINT=(
sin23_0300 sin23_0308 sin23_0316 sin23_0324 sin23_0332 sin23_0340 sin23_0348 sin23_0356 sin23_0364 sin23_0372
sin23_0380 sin23_0388 sin23_0396 sin23_0404 sin23_0412 sin23_0420 sin23_0428 sin23_0436 sin23_0444 sin23_0452
sin23_0460 sin23_0468 sin23_0476 sin23_0484 sin23_0492 sin23_0500 sin23_0508 sin23_0516 sin23_0524 sin23_0532
sin23_0540 sin23_0548 sin23_0556 sin23_0564 sin23_0572 sin23_0580 sin23_0588 sin23_0596 sin23_0604 sin23_0612
sin23_0620 sin23_0628 sin23_0636 sin23_0644 sin23_0652 sin23_0660 sin23_0668 sin23_0676 sin23_0684 sin23_0692
)

P=${POINT[$SLURM_ARRAY_TASK_ID]}

gundamFitter -c config.yaml -of ./inputs/profile/$P.yaml -O /fitterEngineConfig/generateSamplePlots=false --out-dir outputs/profile/ -t 16
