#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Usage: ./run-js-data.sh <phoetmin> <phoetmax> <jetptmin> <trkptmin> <gammaxi> [pp only] [draw pbpb/pp]"
  echo "Example: ./run-js-data.sh 60 1000 30 1 0"
  exit 1
fi

set -x

echo "compiling macros..."
g++ jetshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetshape || exit 1
g++ draw_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_js || exit 1
g++ plot_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_js || exit 1

if [[ $6 -ne 1 ]]; then
    echo running on pbpb data
    ./jetshape /export/d00/scratch/biran/photon-jet-track/PbPb-Data-skim-170405.root pbpbdata 0 20 $1 $2 $3 recoreco $4 $5 &
    ./jetshape /export/d00/scratch/biran/photon-jet-track/PbPb-Data-skim-170405.root pbpbdata 20 60 $1 $2 $3 recoreco $4 $5 &
    ./jetshape /export/d00/scratch/biran/photon-jet-track/PbPb-Data-skim-170405.root pbpbdata 60 100 $1 $2 $3 recoreco $4 $5 &
    ./jetshape /export/d00/scratch/biran/photon-jet-track/PbPb-Data-skim-170405.root pbpbdata 100 200 $1 $2 $3 recoreco $4 $5 &
    wait

    hadd -f pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root pbpbdata_recoreco_${1}_${3}_*_*.root
    rm pbpbdata_recoreco_${1}_${3}_*_*.root

    ./run-js-plot.sh $1 $2 $3 $4 $5 pbpbdata pbpbdata recoreco
fi

echo running on pp data
./jetshape /export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170405.root ppdata 0 20 $1 $2 $3 srecoreco $4 $5 &
./jetshape /export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170405.root ppdata 20 60 $1 $2 $3 srecoreco $4 $5 &
./jetshape /export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170405.root ppdata 60 100 $1 $2 $3 srecoreco $4 $5 &
./jetshape /export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170405.root ppdata 100 200 $1 $2 $3 srecoreco $4 $5 &
wait

hadd -f ppdata_${1}_${3}_gxi${5}_srecoreco_js.root ppdata_srecoreco_${1}_${3}_*_*.root
rm ppdata_srecoreco_${1}_${3}_*_*.root

./run-js-plot.sh $1 $2 $3 $4 $5 ppdata ppdata srecoreco

if [[ $7 -eq 1 ]]; then
    ./run-js-plot.sh $1 $2 $3 $4 $5 data data data
fi
