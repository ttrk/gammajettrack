#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-ff-shape-closure.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-ff-shape-closure.sh 80 1000 40 1 0 1 pbpbmc ffclosure sgengen sgenreco recogen recoreco"
  exit 1
fi

phoetMin=$1
phoetMax=$2
jetptMin=$3
trkptMin=$4
gammaxi=$5
obs=$6
sample=$7
label=$8
recogenTypes=${@:9}

echo "phoetMin = $phoetMin"
echo "phoetMax = $phoetMax"
echo "jetptMin = $jetptMin"
echo "trkptMin = $trkptMin"
echo "gammaxi  = $gammaxi"
echo "obs      = $obs"
echo "sample   = $sample"
echo "label    = $label"
echo "recogenTypes = $recogenTypes"

if [ $7 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
elif [ $7 = "ppmc" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/pp-MC-skim-170911.root"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1

set -x

echo running closure histograms
for recogen in $recogenTypes; do
  if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${recogen}_ffjs.root ]; then
    hiBinMins=(0  20 60  100)
    hiBinMaxs=(20 60 100 200)
    for i1 in ${!hiBinMins[*]}
    do
      hiBinMin=${hiBinMins[i1]}
      hiBinMax=${hiBinMaxs[i1]}
     ./jetffshape.exe $SKIM $7 $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $recogen $trkptMin $gammaxi $label 0 $obs &
    done
  fi
done
wait

for recogen in $recogenTypes; do
  if [ ! -f ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${recogen}_ffjs.root ]; then
    hadd -f ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${recogen}_ffjs.root ${label}_${sample}_${recogen}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
    rm ${label}_${sample}_${recogen}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
  fi
done

#./run-ff-shape-plot.sh $@

outDir="/export/d00/scratch/"$USER"/GJT-out/results/closure/"
mkdir -p $outDir
mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root $outDir
mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root $outDir
mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root $outDir
