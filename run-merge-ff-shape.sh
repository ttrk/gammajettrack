#!/bin/bash

if [ $# -lt 11 ]; then
  echo "Usage: ./run-merge-ff-shape.sh [inputDir] [outDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-merge-ff-shape.sh /mnt/hadoop/cms/store/user/${USER}/GJT-out/results/ /export/d00/scratch/${USER}/GJT-out/results/ 80 1000 40 1 0 1 pbpbmc ffclosure sgengen sgenreco recogen recoreco"
  exit 1
fi

inputDir=$1
outDir=$2
phoetMin=$3
phoetMax=$4
jetptMin=$5
trkptMin=$6
gammaxi=$7
obs=$8
sample=$9
label=${10}
recogenLevels=${@:11}

echo "inputDir = $inputDir"
echo "outDir   = $outDir"
echo "phoetMin = $phoetMin"
echo "phoetMax = $phoetMax"
echo "jetptMin = $jetptMin"
echo "trkptMin = $trkptMin"
echo "gammaxi  = $gammaxi"
echo "obs      = $obs"
echo "sample   = $sample"
echo "label    = $label"
echo "recogenLevels = $recogenLevels"

set -x

mkdir -p $outDir
for rgLevel in $recogenLevels; do
  if [ ! -f ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
    hadd -f ${outDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ${inputDir}/${label}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
    # rm ${label}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
  fi
done

#./run-ff-shape-plot.sh $@

#outDir="/export/d00/scratch/"$USER"/GJT-out/results/closure/"
#mkdir -p $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root $outDir
