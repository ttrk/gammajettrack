#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: ./save-dphidetarefrecoJet-ptBin.sh [inputFile] [outputFile] [sample]"
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
