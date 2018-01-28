#!/bin/bash

if [ $# -lt 8 ]; then
  echo "Usage: ./run-js-closure.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
  echo "Example: ./run-js-closure.sh 60 1000 30 1 0 pbpbmc closure sgengen sgenreco recogen recoreco"
  exit 1
fi

ARGS=()

while [ $# -gt 0 ]; do
    case "$1" in
        -g)             GROUP="$2"; shift 2 ;;
        --group=*)      GROUP="${1#*=}"; shift ;;
        -j)             JOBS="$2"; shift 2 ;;
        --jobs=*)       JOBS="${1#*=}"; shift ;;
        -n)             NICE="$2"; shift 2 ;;
        --nice=*)       NICE="${1#*=}"; shift ;;
        -*|--*)         [[ $1 =~ ^-?[0-9]+$ ]] && \
                            { ARGS+=("$1"); shift; } || \
                            { echo -e "invalid option: $1\n"; exit 1; } ;;
        *)              ARGS+=("$1"); shift ;;
    esac
done

set -- "${ARGS[@]}"

if [ $6 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-180115.root"
    TOTAL=51
elif [ $6 = "ppmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-180115.root"
    TOTAL=15
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
make jetshape

GROUP=${GROUP:-"def"}
JOBS=${JOBS:-"+0"}

FLAGS="--linebuffer"

[ -n "$NICE" ] && PREFIX+="nice -n $NICE "

echo running closure histograms
for i in ${@:8}; do
    for slice in $(seq 0 $TOTAL); do
        $PREFIX sem --id rjc-$GROUP -j$JOBS $FLAGS  \
                "./jetshape $SKIM $6                \
                0 20 $1 $2 $3 $i $4 $5 $7 0 $slice"
        $PREFIX sem --id rjc-$GROUP -j$JOBS $FLAGS  \
                "./jetshape $SKIM $6                \
                20 60 $1 $2 $3 $i $4 $5 $7 0 $slice"
        $PREFIX sem --id rjc-$GROUP -j$JOBS $FLAGS  \
                "./jetshape $SKIM $6                \
                60 100 $1 $2 $3 $i $4 $5 $7 0 $slice"
        $PREFIX sem --id rjc-$GROUP -j$JOBS $FLAGS  \
                "./jetshape $SKIM $6                \
                100 200 $1 $2 $3 $i $4 $5 $7 0 $slice"
    done
done
sem --id rjc-$GROUP --wait

for i in ${@:8}; do
    hadd -f ${7}_${6}_${1}_${3}_gxi${5}_${i}_js.root ${7}_${6}_${i}_${1}_${3}_${5}_*_*_*.root
    rm ${7}_${6}_${i}_${1}_${3}_${5}_*_*.root
done

./run-js-arithmetic.sh $@
