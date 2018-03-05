#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-ff-shape-data.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-ff-shape-data.sh /export/d00/scratch/$USER/GJT-out/results/closure/ 80 1000 40 1 0 1 pbpbdata jsdata sgengen sgenreco recogen recoreco"
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

if [ $sample = "pbpbdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
elif [ $sample = "ppdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/pp-Data-skim-170911.root"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1

mkdir -p $outputDir

set -x

echo running closure histograms
for rgLevel in $recogenLevels; do
#  if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
    hiBinMins=(0  20 60  100 0  60)
    hiBinMaxs=(20 60 100 200 60 200)
    if [ $sample = "ppdata" ] && [ $rgLevel = "recoreco" ]; then
      hiBinMins=(100)
      hiBinMaxs=(200)
    fi
    for i1 in ${!hiBinMins[*]}
    do
      hiBinMin=${hiBinMins[i1]}
      hiBinMax=${hiBinMaxs[i1]}
     ./jetffshape.exe $SKIM $sample $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $rgLevel $trkptMin $gammaxi ${outputDir}/$label 0 $obs &
    done
#  fi
done
wait

for rgLevel in $recogenLevels; do
  outputFile=${outputDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root
  inputPrefix=${outputDir}/${label}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}
#  if [ ! -f $outputFile ]; then
  hadd -f $outputFile ${inputPrefix}_*_*.root
  rm ${inputPrefix}_*_*.root
#  fi
done

outputFilesMerged=${outputDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
hadd -f $outputFilesMerged ${outputDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root

#./run-ff-shape-plot.sh $@
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root $outDir
