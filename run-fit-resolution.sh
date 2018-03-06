#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-fit-resolution.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-fit-resolution.sh /export/d00/scratch/"$USER"/GJT-out/results/ 60 1000 30 1 0 2 pbpbmc jsclosure reco0gen0 sref0gen0 sptref0gen0 sphiref0gen0 setaref0gen0"
  exit 1
fi

outputDir=$1
phoetMin=$2
phoetMax=$3
jetptMin=$4
trkptMin=$5
gammaxi=$6
obs=$7
sample=$8
label=$9
recogenLevels=${@:10}

echo "outputDir = $outputDir"
echo "phoetMin  = $phoetMin"
echo "phoetMax  = $phoetMax"
echo "jetptMin  = $jetptMin"
echo "trkptMin  = $trkptMin"
echo "gammaxi   = $gammaxi"
echo "obs       = $obs"
echo "sample    = $sample"
echo "label     = $label"
echo "recogenLevels = $recogenLevels"

mkdir -p $outputDir

g++ fit_resolution.C $(root-config --cflags --libs) -Werror -Wall -O2 -o fit_resolution.exe || exit 1
echo "g++ fit_resolution.C $(root-config --cflags --libs) -Werror -Wall -O2 -o fit_resolution.exe || exit 1"

if [[ ${label} == "jsclosure" ]]; then
    echo "running ${label}"
    set -x
    inputFile=${outputDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
    if [ ! -f $inputFile ]; then
      echo "inputFile : $inputFile does not exist. Exiting."
      exit 1
    fi
    outputFile="${inputFile/_merged.root/_fit_res.root}"

    ./fit_resolution.exe $sample $inputFile $outputFile $phoetMin $recogenLevels
fi
