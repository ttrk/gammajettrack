#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-ff-closure.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [defnFF] [sample] [label] [types...]"
  echo "Example: ./run-ff-closure.sh 80 1000 40 1 0 1 pbpbmc ffclosure sgengen sgenreco recogen recoreco"
  exit 1
fi

echo "phoetmin = $1"
echo "phoetmax = $2"
echo "jetptmin = $3"
echo "trkptmin = $4"
echo "gammaxi  = $5"
echo "defnFF   = $6"
echo "sample   = $7"
echo "label    = $8"
echo "types    = $9"

if [ $7 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170824.root"
elif [ $7 = "ppmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-170824.root"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetff.exe || exit 1

set -x

echo running closure histograms
for i in ${@:9}; do
  if [ ! -f ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_${i}_ff.root ]; then
    ./jetff.exe $SKIM $7 0 20 $1 $2 $3 $i $4 $5 $8 0 $6 &
    ./jetff.exe $SKIM $7 20 60 $1 $2 $3 $i $4 $5 $8 0 $6 &
    ./jetff.exe $SKIM $7 60 100 $1 $2 $3 $i $4 $5 $8 0 $6 &
    ./jetff.exe $SKIM $7 100 200 $1 $2 $3 $i $4 $5 $8 0 $6 &
  fi
done
wait

for i in ${@:9}; do
  if [ ! -f ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_${i}_ff.root ]; then
    hadd -f ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_${i}_ff.root ${8}_${7}_${i}_${1}_${3}_${5}_${6}_*_*.root
    rm ${8}_${7}_${i}_${1}_${3}_${5}_${6}_*_*.root
  fi
done

./run-ff-plot.sh $@

outDir="/export/d00/scratch/"$USER"/GJT-out/results/closure/"
mkdir -p $outDir
mv ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_*_ff.root $outDir
mv ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root $outDir
mv ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root $outDir
