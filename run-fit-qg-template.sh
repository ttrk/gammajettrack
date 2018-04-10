#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./run-fit-qg-template.sh [inputMC] [inputData] [output]"
  echo "Example: ./run-fit-qg-template.sh /export/d00/scratch/tatar/GJT-out/results/closureQG/jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_merged.root /export/d00/scratch/tatar/GJT-out/results/jsdata_pbpbdata_60_30_gxi0_obs2_ffjs_final.root ./fit_qg_template_pbpb_gxi0_obs2.root"
  exit 1
fi

inputMC=$1
inputData=$2
output=$3

echo "inputMC = $inputMC"
echo "inputData = $inputData"
echo "output = $output"

g++ fit_qg_template.C $(root-config --cflags --libs) -Werror -Wall -O2 -o fit_qg_template.exe || exit 1

set -x

./fit_qg_template.exe $inputMC $inputData $output

