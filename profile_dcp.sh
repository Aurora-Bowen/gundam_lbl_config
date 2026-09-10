#!/bin/bash

#SBATCH --job-name=HK_prof_dcp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=0-47%8
#SBATCH --output=outputs/logs/profile_dcp_%A_%a.out
#SBATCH --mem-per-cpu=4G

# Profile likelihood in delta-CP. Forty-eight independent fits, one per array
# task, each freezing delta-CP at one grid point and minimising the other five
# oscillation parameters and all 70 systematics around it.
#
# Grid: -180 to +172.5 degrees in steps of 7.5, over the whole cyclicRange
# rather than a few sigma either side, since delta-CP is cyclic and its interval
# is quoted over all of it. Prior -1.601 rad, which is -91.7 degrees.
#
# 48 points and not 50: 7.5 degrees puts 0 and both maximal-CP-violation points
# on the grid, where 360/50 = 7.2 would miss +-90. Forty-eight steps of 7.5
# close 360 exactly, so +180 is not a separate task - it is task 0, at -180.
#
# Every task writes its own file into outputs/profile/, named from the override
# it was given, so they cannot collide. See inputs/profile/README.md.
#
#   sbatch profile_dcp.sh        from this directory

source "$HOME/setup_neut.sh"

[ -n "$data_dir" ] || { echo "data_dir is not set"; exit 1; }

mkdir -p outputs/profile/

POINT=(
dcp_m1800 dcp_m1725 dcp_m1650 dcp_m1575 dcp_m1500 dcp_m1425 dcp_m1350 dcp_m1275
dcp_m1200 dcp_m1125 dcp_m1050 dcp_m0975 dcp_m0900 dcp_m0825 dcp_m0750 dcp_m0675
dcp_m0600 dcp_m0525 dcp_m0450 dcp_m0375 dcp_m0300 dcp_m0225 dcp_m0150 dcp_m0075
dcp_0000  dcp_p0075 dcp_p0150 dcp_p0225 dcp_p0300 dcp_p0375 dcp_p0450 dcp_p0525
dcp_p0600 dcp_p0675 dcp_p0750 dcp_p0825 dcp_p0900 dcp_p0975 dcp_p1050 dcp_p1125
dcp_p1200 dcp_p1275 dcp_p1350 dcp_p1425 dcp_p1500 dcp_p1575 dcp_p1650 dcp_p1725
)

P=${POINT[$SLURM_ARRAY_TASK_ID]}

gundamFitter -c config.yaml -of ./inputs/profile/$P.yaml -O /fitterEngineConfig/generateSamplePlots=false --out-dir outputs/profile/ -t 8
