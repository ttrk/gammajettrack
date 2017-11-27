#!/bin/bash

if [ $# -ne 7 ]; then
  echo "Usage: ./run-ff-shape-data.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [defnFF] [coll]"
  echo "       [coll] - pbpb, pp, both"
  echo "Example: ./run-ff-shape-data.sh 80 1000 40 1 0 0 2"
  exit 1
fi

echo "phoetmin = $1"
echo "phoetmax = $2"
echo "jetptmin = $3"
echo "trkptmin = $4"
echo "gammaxi  = $5"
echo "defnFF   = $6"
echo "coll     = $7"

g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1
echo "g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1"

PBPBSKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
PPSKIM="/export/d00/scratch/tatar/GJT-out/pp-Data-skim-170911.root"

outDir="/export/d00/scratch/"$USER"/GJT-out/results/"
mkdir -p $outDir
echo "outDir : $outDir"

systematics=13
set +x

if [[ $7 == "pbpb" ]]; then
    echo "running on pbpb data"
    echo "PBPBSKIM : $PBPBSKIM"
    set -x
    outPrefix=$outDir"data"
    ./jetffshape.exe $PBPBSKIM pbpbdata 0 20 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PBPBSKIM pbpbdata 20 60 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PBPBSKIM pbpbdata 60 100 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PBPBSKIM pbpbdata 100 200 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PBPBSKIM pbpbdata 0 60 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PBPBSKIM pbpbdata 60 200 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    wait

    hadd -f ${outPrefix}_pbpbdata_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root ${outPrefix}_pbpbdata_recoreco_${1}_${3}_${5}_${6}_*_*.root
    rm ${outPrefix}_pbpbdata_recoreco_${1}_${3}_${5}_${6}_*_*.root

    ./run-ff-shape-plot.sh $1 $2 $3 $4 $5 $6 pbpbdata data recoreco
fi

if [[ $7 == "pp" ]]; then
    echo "running on pp data"
    echo "PPSKIM : $PPSKIM"
    set -x
    outPrefix=$outDir"data"
    ./jetffshape.exe $PPSKIM ppdata 0 20 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 20 60 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 60 100 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 100 200 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 0 60 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 60 200 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &

    ./jetffshape.exe $PPSKIM ppdatareweight 0 20 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdatareweight 20 60 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdatareweight 60 100 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdatareweight 100 200 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdatareweight 0 60 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdatareweight 60 200 $1 $2 $3 srecoreco $4 $5 $outPrefix $systematics $6 &

    ./jetffshape.exe $PPSKIM ppdata 0 20 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 20 60 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 60 100 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 100 200 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 0 60 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    ./jetffshape.exe $PPSKIM ppdata 60 200 $1 $2 $3 recoreco $4 $5 $outPrefix $systematics $6 &
    wait

    hadd -f ${outPrefix}_ppdata_${1}_${3}_gxi${5}_defnFF${6}_srecoreco_ff.root ${outPrefix}_ppdata_srecoreco_${1}_${3}_${5}_${6}_*_*.root
    hadd -f ${outPrefix}_ppdatareweight_${1}_${3}_gxi${5}_defnFF${6}_srecoreco_ff.root ${outPrefix}_ppdatareweight_srecoreco_${1}_${3}_${5}_${6}_*_*.root
    hadd -f ${outPrefix}_ppdata_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root ${outPrefix}_ppdata_recoreco_${1}_${3}_${5}_${6}_*_*.root
    rm ${outPrefix}_ppdata_srecoreco_${1}_${3}_${5}_${6}_*_*.root
    rm ${outPrefix}_ppdatareweight_srecoreco_${1}_${3}_${5}_${6}_*_*.root
    rm ${outPrefix}_ppdata_recoreco_${1}_${3}_${5}_${6}_*_*.root

    ./run-ff-shape-plot.sh $1 $2 $3 $4 $5 $6 ppdata data 0
fi

set +x
if [[ $7 == "both" ]]; then
    echo "running plotting for pbpb and pp"
    set -x
    ./run-ff-shape-plot.sh $1 $2 $3 $4 $5 $6 data data data
fi
