#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./save-dphidetarefrecoJet-ptBin.sh [inputFile] [outputFile] [sample]"
  echo "Example: ./run-save-dphidetarefrecoJet-ptBin.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_ppmc_60_30_gxi0_obs2_ffjs_merged.root pp-weights.root ppmc"
  echo "Example: ./run-save-dphidetarefrecoJet-ptBin.sh /export/d00/scratch/${USER}/GJT-out/results/closure/jsclosure_pbpbmc_60_30_gxi0_obs2_ffjs_merged.root PbPb-weights.root pbpbmc"
  exit 1
fi

inputFile=$1
outputFile=$2
sample=$3

echo "inputFile  = $inputFile"
echo "outputFile = $outputFile"
echo "sample     = $sample"

g++ save_dphidetarefrecoJet_ptBin.C $(root-config --cflags --libs) -Werror -Wall -O2 -o save_dphidetarefrecoJet_ptBin.exe || exit 1
echo "g++ save_dphidetarefrecoJet_ptBin.C $(root-config --cflags --libs) -Werror -Wall -O2 -o save_dphidetarefrecoJet_ptBin.exe || exit 1"

./save_dphidetarefrecoJet_ptBin.exe $inputFile $outputFile $sample
