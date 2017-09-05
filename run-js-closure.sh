#!/bin/bash

if [ $# -lt 8 ]; then
  echo "Usage: ./run-js-closure.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
  echo "Example: ./run-js-closure.sh 80 1000 40 1 0 pbpbmc jsclosure sgengen sgenreco recogen recoreco"
  exit 1
fi

if [ $6 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170824.root"
elif [ $6 = "ppmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-170824.root"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetshape || exit 1

set -x

echo running closure histograms
for i in ${@:8}; do
  if [ ! -f ${7}_${6}_${1}_${3}_gxi${5}_${i}_js.root ]; then
    ./jetshape $SKIM $6 0 20 $1 $2 $3 $i $4 $5 $7 &
    ./jetshape $SKIM $6 20 60 $1 $2 $3 $i $4 $5 $7 &
    ./jetshape $SKIM $6 60 100 $1 $2 $3 $i $4 $5 $7 &
    ./jetshape $SKIM $6 100 200 $1 $2 $3 $i $4 $5 $7 &
  fi
done
wait

for i in ${@:8}; do
  if [ ! -f ${7}_${6}_${1}_${3}_gxi${5}_${i}_js.root ]; then
    hadd -f ${7}_${6}_${1}_${3}_gxi${5}_${i}_js.root ${7}_${6}_${i}_${1}_${3}_${5}_*_*.root
    rm ${7}_${6}_${i}_${1}_${3}_${5}_*_*.root
  fi
done

./run-js-plot.sh $@
