#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./run-plot-js-corrections.sh [inputFile] [outputDir] [sample]"
  echo "Example: ./run-plot-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_ppmc_60_30_gxi0_obs2_ffjs_merged.root jscorrections_ppmc_60_30_gxi0_obs2.root ppmc"
  echo "Example: ./run-plot-js-corrections.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_merged.root  jscorrections_pbpbmc_60_30_gxi0_obs2.root pbpbmc"
  exit 1
fi

inputFile=$1
outputDir=$2
sample=$3

echo "inputFile = $inputFile"
echo "outputDir = $outputDir"
echo "sample    = $sample"

g++ plot_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_js_corrections.exe || exit 1
echo "g++ plot_js_corrections.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_js_corrections.exe || exit 1"

mkdir -p $outputDir
./plot_js_corrections.exe $inputFile $outputDir $sample
