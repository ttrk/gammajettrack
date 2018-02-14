#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./run-calc-js-corrections.sh [inputFile] [outputFile] [sample]"
  echo "Example: ./run-calc-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_ppmc_60_30_gxi0_obs2_ffjs_merged.root jscorrections_ppmc_60_30_gxi0_obs2.root ppmc"
  echo "Example: ./run-calc-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_merged.root  jscorrections_pbpbmc_60_30_gxi0_obs2.root pbpbmc"
  exit 1
fi

inputFile=$1
outputFile=$2
sample=$3

echo "inputFile  = $inputFile"
echo "outputFile = $outputFile"
echo "sample     = $sample"

g++ calc_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_js_corrections.exe || exit 1
echo "g++ calc_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_js_corrections.exe || exit 1"

./calc_js_corrections.exe $inputFile $outputFile $sample
