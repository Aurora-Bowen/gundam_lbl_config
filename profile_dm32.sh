#!/bin/bash

#SBATCH --job-name=HK_prof_dm32
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-49%8
#SBATCH --output=outputs/logs/profile_dm32_%A_%a.out
#SBATCH --mem-per-cpu=4G

# Profile likelihood in dm^2_32. Fifty independent fits, one per array task,
# each freezing dm^2_32 at one grid point and minimising the other five
# oscillation parameters and all 70 systematics around it.
#
# Grid: 2.000E-3 to 2.980E-3 in steps of 2.0E-5, which is 0.59 sigma per step
# against the 3.4E-5 in MakeCovariances.C. Prior 2.509E-3, upper parameterLimit
# 3.000E-3.
#
# Every task writes its own file into outputs/profile/, named from the override
# it was given, so they cannot collide. See inputs/profile/README.md.
#
#   sbatch profile_dm32.sh       from this directory

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/profile/

POINT=(
dm32_2000 dm32_2020 dm32_2040 dm32_2060 dm32_2080 dm32_2100 dm32_2120 dm32_2140 dm32_2160 dm32_2180
dm32_2200 dm32_2220 dm32_2240 dm32_2260 dm32_2280 dm32_2300 dm32_2320 dm32_2340 dm32_2360 dm32_2380
dm32_2400 dm32_2420 dm32_2440 dm32_2460 dm32_2480 dm32_2500 dm32_2520 dm32_2540 dm32_2560 dm32_2580
dm32_2600 dm32_2620 dm32_2640 dm32_2660 dm32_2680 dm32_2700 dm32_2720 dm32_2740 dm32_2760 dm32_2780
dm32_2800 dm32_2820 dm32_2840 dm32_2860 dm32_2880 dm32_2900 dm32_2920 dm32_2940 dm32_2960 dm32_2980
)

P=${POINT[$SLURM_ARRAY_TASK_ID]}

gundamFitter -c config.yaml -of ./inputs/profile/$P.yaml -O /fitterEngineConfig/generateSamplePlots=false --out-dir outputs/profile/ -t 16
