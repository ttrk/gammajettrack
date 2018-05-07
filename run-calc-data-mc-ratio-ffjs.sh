#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./run-calc-data-mc-ratio-ffjs.sh [inputData] [inputMC] [output]"
  echo "Example: ./run-calc-data-mc-ratio-ffjs.sh /export/d00/scratch/tatar/GJT-out/results//jsdata_ppdata_60_30_gxi0_obs2_ffjs_final.root /export/d00/scratch/tatar/GJT-out/results/closure/jsclosure_ppmc_60_30_gxi0_obs2_ffjs_final.root /export/d00/scratch/tatar/GJT-out/results/jsdata_pp_data_mc_ratio_60_30_gxi0_obs2.root"
  exit 1
fi

inputData=$1
inputMC=$2
output=$3

echo "inputData = $inputData"
echo "inputMC = $inputMC"
echo "output = $output"

g++ calc_data_mc_ratio.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_data_mc_ratio.exe || exit 1

set -x

./calc_data_mc_ratio.exe $inputData $inputMC $output

