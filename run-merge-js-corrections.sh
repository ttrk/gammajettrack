#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./run-merge-js-corrections.sh [inputFile] [outputFile] [sample]"
  echo "Example: ./run-merge-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_ppmc_60_30_gxi0_obs2_ffjs_merged.root merged_ppmc_60_30_gxi0_obs2.root ppmc"
  echo "Example: ./run-merge-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_merged.root merged_pbpbmc_60_30_gxi0_obs2.root pbpbmc"
  exit 1
fi

inputFile=$1
outputFile=$2
sample=$3

echo "inputFile  = $inputFile"
echo "outputFile = $outputFile"
echo "sample     = $sample"

g++ merge_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o merge_js_corrections.exe || exit 1
echo "g++ merge_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o merge_js_corrections.exe || exit 1"

./merge_js_corrections.exe $inputFile $outputFile $sample
