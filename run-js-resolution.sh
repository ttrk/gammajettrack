#!/bin/bash

if [ $# -lt 7 ]; then
  echo "Usage: ./run-js-resolution.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label]"
  echo "Example: ./run-js-resolution.sh 60 1000 30 1 0 pbpbmc resolution"
  exit 1
fi

if [ $6 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
    BKGSKIM=""
elif [ $6 = "ppmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-180115.root"
    BKGSKIM=""
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
make jetres || exit 1

set -x

echo running resolution histograms
./jetres $SKIM "$BKGSKIM" $6 0 20 $1 $2 $3 a $4 $5 $7 &
./jetres $SKIM "$BKGSKIM" $6 20 60 $1 $2 $3 a $4 $5 $7 &
./jetres $SKIM "$BKGSKIM" $6 60 100 $1 $2 $3 a $4 $5 $7 &
./jetres $SKIM "$BKGSKIM" $6 100 200 $1 $2 $3 a $4 $5 $7 &
wait
