#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "Usage: ./plot-js-systematics.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi]"
    exit 1
fi

# compile plot_js

SYSTEMATICS=(jes_up_plus_plus_plus jer_minus_minus_plus_plus pes purity_up_plus_plus_plus tracking iso)
SYSLABEL=(JES JER "photon energy" "photon purity" tracking "photon isolation")
SAMPLE=(pbpb pp)
TYPE=(recoreco srecoreco)
SMPLLABEL=(PbPb pp)

g++ plot_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_js || exit 1

set -x

for SMPL in 0 1
do
    for SYS in {0..5}
    do
        HISTLIST=${SYSTEMATICS[SYS]}_${1}_${2}_${3}_gxi${5}.list
        if [ -f $HISTLIST ]; then
            rm $HISTLIST
        fi
        touch $HISTLIST

        echo -e "${SMPLLABEL[SMPL]} data" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_0_20_${SYSTEMATICS[SYS]}_nominal" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_20_60_${SYSTEMATICS[SYS]}_nominal" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_60_100_${SYSTEMATICS[SYS]}_nominal" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_100_200_${SYSTEMATICS[SYS]}_nominal" >> $HISTLIST
        echo -e "${SYSLABEL[SYS]}" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_0_20_${SYSTEMATICS[SYS]}_variation" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_20_60_${SYSTEMATICS[SYS]}_variation" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_60_100_${SYSTEMATICS[SYS]}_variation" >> $HISTLIST
        echo -e "hjs_final_${SAMPLE[SMPL]}data_${TYPE[SMPL]}_100_200_${SYSTEMATICS[SYS]}_variation" >> $HISTLIST

        ./plot_js /export/d00/scratch/biran/photon-jet-track/js_systematics/data_${1}_${3}_gxi${5}-systematics.root ${SYSTEMATICS[SYS]}_${SAMPLE[SMPL]}_plot_${1}_${2}_${3}_gxi${5} $HISTLIST 1 $5 $1 $3

        rm $HISTLIST
    done
done
