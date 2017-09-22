#!/bin/bash

if [ $# -ne 7 ]; then
  echo -e "Usage: ./run-js-data.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label]"
  echo -e "Example: ./run-js-data.sh 60 1000 30 1 0 both nominal"
  exit 1
fi

echo -e "compiling macros..."
make jetshape

PBPBSKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
PPSKIM="/export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170911.root"

set -x

if [[ $6 != "pp" ]]; then
    echo -e "running on pbpb data"
    ./jetshape $PBPBSKIM pbpbdata 0 20 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PBPBSKIM pbpbdata 20 60 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PBPBSKIM pbpbdata 60 100 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PBPBSKIM pbpbdata 100 200 $1 $2 $3 recoreco $4 $5 $7 &
    wait

    hadd -f ${7}_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root ${7}_pbpbdata_recoreco_${1}_${3}_${5}_*_*.root
    rm ${7}_pbpbdata_recoreco_${1}_${3}_${5}_*_*.root

    ./run-js-arithmetic.sh ${@:1:5} pbpbdata ${7} recoreco
fi

if [[ $6 != "pbpb" ]]; then
    echo -e "running on pp data"
    ./jetshape $PPSKIM ppdata 0 20 $1 $2 $3 srecoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 20 60 $1 $2 $3 srecoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 60 100 $1 $2 $3 srecoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 100 200 $1 $2 $3 srecoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 0 20 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 20 60 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 60 100 $1 $2 $3 recoreco $4 $5 $7 &
    ./jetshape $PPSKIM ppdata 100 200 $1 $2 $3 recoreco $4 $5 $7 &
    wait

    hadd -f ${7}_ppdata_${1}_${3}_gxi${5}_srecoreco_js.root ${7}_ppdata_srecoreco_${1}_${3}_${5}_*_*.root
    hadd -f ${7}_ppdata_${1}_${3}_gxi${5}_recoreco_js.root ${7}_ppdata_recoreco_${1}_${3}_${5}_*_*.root
    rm ${7}_ppdata_srecoreco_${1}_${3}_${5}_*_*.root
    rm ${7}_ppdata_recoreco_${1}_${3}_${5}_*_*.root

    ./run-js-arithmetic.sh ${@:1:5} ppdata $7 srecoreco
    ./run-js-arithmetic.sh ${@:1:5} ppdata $7 recoreco

    ./run-js-arithmetic.sh ${@:1:5} ppdata $7 recoreco srecoreco
fi

if [[ $6 == "both" ]]; then
    ./run-js-arithmetic.sh ${@:1:5} ${7} ${7} data
fi
