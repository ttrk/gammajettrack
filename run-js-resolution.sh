#!/bin/bash

helpmsg() {
    echo -e 'usage:   ./run-js-resolution.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label]'
    echo -e 'example: ./run-js-resolution.sh 60 1000 30 1 0 pbpbmc resolution\n'
    echo -e '   -f, --fit       refit only'
    echo -e '   -g, --group     group identifier'
    echo -e '   -h, --help      show (this) help message'
    echo -e '   -j, --jobs      jobs relative to number of cores'
    echo -e '   -n, --nice      niceness\n'
}

ARGS=()

while [ $# -gt 0 ]; do
    case "$1" in
        -f|--fit)       fitonly=1; shift ;;
        -g)             GROUP="$2"; shift 2 ;;
        --group=*)      GROUP="${1#*=}"; shift ;;
        -h|--help)      helpmsg; exit ;;
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

[ $# -lt 7 ] && { helpmsg; exit; }

case "$6" in
    pbpbmc)
        SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-180115.root"
        TOTAL=51
        ;;
    ppmc)
        SKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-180115.root"
        TOTAL=15
        ;;
    *)
        echo "invalid sample"
        exit 1
        ;;
esac

echo "compiling macros..."
make jetres fitjetres || exit 1

GROUP=${GROUP:-"def"}
JOBS=${JOBS:-"+0"}

FLAGS="--linebuffer"

[ -n "$NICE" ] && PREFIX+="nice -n $NICE "

if [ ! -n "$fitonly" ]; then
  echo running resolution histograms
    for slice in $(seq 0 $TOTAL); do
        $PREFIX sem --id rjr-$GROUP -j$JOBS $FLAGS  \
                "./jetres $SKIM $6 0 20 $1 $2 $3 a $4 $5 $7 0 $slice"
        $PREFIX sem --id rjr-$GROUP -j$JOBS $FLAGS  \
                "./jetres $SKIM $6 20 60 $1 $2 $3 a $4 $5 $7 0 $slice"
        $PREFIX sem --id rjr-$GROUP -j$JOBS $FLAGS  \
                "./jetres $SKIM $6 60 100 $1 $2 $3 a $4 $5 $7 0 $slice"
        $PREFIX sem --id rjr-$GROUP -j$JOBS $FLAGS  \
                "./jetres $SKIM $6 100 200 $1 $2 $3 a $4 $5 $7 0 $slice"
    done
    sem --id rjr-$GROUP --wait

    hadd -f ${7}_${6}_${1}_${3}_${5}_0_20.root ${7}_${6}_${1}_${3}_${5}_0_20_*.root
    hadd -f ${7}_${6}_${1}_${3}_${5}_20_60.root ${7}_${6}_${1}_${3}_${5}_20_60_*.root
    hadd -f ${7}_${6}_${1}_${3}_${5}_60_100.root ${7}_${6}_${1}_${3}_${5}_60_100_*.root
    hadd -f ${7}_${6}_${1}_${3}_${5}_100_200.root ${7}_${6}_${1}_${3}_${5}_100_200_*.root

    rm ${7}_${6}_${1}_${3}_${5}_0_20_*.root
    rm ${7}_${6}_${1}_${3}_${5}_20_60_*.root
    rm ${7}_${6}_${1}_${3}_${5}_60_100_*.root
    rm ${7}_${6}_${1}_${3}_${5}_100_200_*.root
fi

./fitjetres ${7}_${6}_${1}_${3}_${5}_0_20.root $6 0 20 ${@:8}
./fitjetres ${7}_${6}_${1}_${3}_${5}_20_60.root $6 20 60 ${@:8}
./fitjetres ${7}_${6}_${1}_${3}_${5}_60_100.root $6 60 100 ${@:8}
./fitjetres ${7}_${6}_${1}_${3}_${5}_100_200.root $6 100 200 ${@:8}
