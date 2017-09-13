#!/bin/bash

if [ $# -ne 6 ]; then
  echo "Usage: ./run-js-data.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [which]"
  echo "       [which] - 0: both, 1: pp, 2: pbpb"
  echo "Example: ./run-js-data.sh 80 1000 40 1 0 0"
  exit 1
fi

echo "compiling macros..."
# g++ jetshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetshape || exit 1
# g++ draw_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_js || exit 1
# g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results || exit 1

PBPBSKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
PPSKIM="/export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170911.root"

set -x

if [[ $6 -ne 1 ]]; then
    echo running on pbpb data
    ./jetshape $PBPBSKIM pbpbdata 0 20 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PBPBSKIM pbpbdata 20 60 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PBPBSKIM pbpbdata 60 100 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PBPBSKIM pbpbdata 100 200 $1 $2 $3 recoreco $4 $5 data &
    wait

    hadd -f data_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root data_pbpbdata_recoreco_${1}_${3}_${5}_*_*.root
    rm data_pbpbdata_recoreco_${1}_${3}_${5}_*_*.root

    ./run-js-plot.sh $1 $2 $3 $4 $5 pbpbdata data recoreco
fi

if [[ $6 -ne 2 ]]; then
    echo running on pp data
    ./jetshape $PPSKIM ppdata 0 20 $1 $2 $3 srecoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 20 60 $1 $2 $3 srecoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 60 100 $1 $2 $3 srecoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 100 200 $1 $2 $3 srecoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 0 20 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 20 60 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 60 100 $1 $2 $3 recoreco $4 $5 data &
    ./jetshape $PPSKIM ppdata 100 200 $1 $2 $3 recoreco $4 $5 data &
    wait

    hadd -f data_ppdata_${1}_${3}_gxi${5}_srecoreco_js.root data_ppdata_srecoreco_${1}_${3}_${5}_*_*.root
    hadd -f data_ppdata_${1}_${3}_gxi${5}_recoreco_js.root data_ppdata_recoreco_${1}_${3}_${5}_*_*.root
    rm data_ppdata_srecoreco_${1}_${3}_${5}_*_*.root
    rm data_ppdata_recoreco_${1}_${3}_${5}_*_*.root

    ./run-js-plot.sh $1 $2 $3 $4 $5 ppdata data srecoreco
    ./run-js-plot.sh $1 $2 $3 $4 $5 ppdata data recoreco
fi

if [[ $6 -eq 0 ]]; then
    ./run-js-plot.sh $1 $2 $3 $4 $5 data data data
fi
