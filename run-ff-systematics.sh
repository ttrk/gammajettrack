#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "Usage: ./run-ff-systematics.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi]"
    exit 1
fi

g++ draw_ff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff.exe || exit 1
g++ calc_ff_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ff_systematics.exe || exit 1

set -x

SAMPLE=(pbpb pp)
SYSTEMATIC=(nominal jer03 jer06 jer10 jes03 jes06 jes10 jesm03 jesm06 jesm10 pes03 pes06 pes10 pesm03 pesm06 pesm10 trk purity_up purity_up_minus purity_down_plus purity_down)

for SMPL in 0 1
do
    for SYS in {0..16}
    do
        ./draw_ff.exe ${SAMPLE[SMPL]}data /export/d00/scratch/biran/photon-jet-track/v3syst/${SYSTEMATIC[SYS]}/all_${1}_${2}_${3}_gammaxi${5}.root /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/${SYSTEMATIC[SYS]}_${1}_${2}_${3}_gammaxi${5}.root ${1} 0 recoreco
    done

    ./draw_ff.exe ${SAMPLE[SMPL]}data /export/d00/scratch/biran/photon-jet-track/v3syst/nominal/all_${1}_${2}_${3}_gammaxi${5}.root /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/purity_up_${1}_${2}_${3}_gammaxi${5}.root ${1} 2 recoreco
    ./draw_ff.exe ${SAMPLE[SMPL]}data /export/d00/scratch/biran/photon-jet-track/v3syst/nominal/all_${1}_${2}_${3}_gammaxi${5}.root /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/purity_up_minus_${1}_${2}_${3}_gammaxi${5}.root ${1} 1 recoreco
    ./draw_ff.exe ${SAMPLE[SMPL]}data /export/d00/scratch/biran/photon-jet-track/v3syst/nominal/all_${1}_${2}_${3}_gammaxi${5}.root /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/purity_down_plus_${1}_${2}_${3}_gammaxi${5}.root ${1} -1 recoreco
    ./draw_ff.exe ${SAMPLE[SMPL]}data /export/d00/scratch/biran/photon-jet-track/v3syst/nominal/all_${1}_${2}_${3}_gammaxi${5}.root /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/purity_down_${1}_${2}_${3}_gammaxi${5}.root ${1} -2 recoreco

    SYSLIST=ff_systematics_${1}_${3}_${5}_${SAMPLE[SMPL]}.list
    if [ -f $SYSLIST ]; then
        rm $SYSLIST
    fi
    touch $SYSLIST

    for SYS in {1..20}
    do
        echo -e "/export/d00/scratch/biran/photon-jet-track/ff_syst_v3/${SYSTEMATIC[SYS]}_${1}_${2}_${3}_gammaxi${5}.root" >> $SYSLIST
    done

    HISTLIST=ff_hist_${1}_${3}_${5}_${SAMPLE[SMPL]}.list
    if [ -f $HISTLIST ]; then
        rm $HISTLIST
    fi
    touch $HISTLIST
    echo -e "hff_final_${SAMPLE[SMPL]}data_recoreco_0_20" >> $HISTLIST
    echo -e "hff_final_${SAMPLE[SMPL]}data_recoreco_20_60" >> $HISTLIST
    echo -e "hff_final_${SAMPLE[SMPL]}data_recoreco_60_100" >> $HISTLIST
    echo -e "hff_final_${SAMPLE[SMPL]}data_recoreco_100_200" >> $HISTLIST

    ./calc_ff_systematics.exe /export/d00/scratch/biran/photon-jet-track/ff_syst_v3/nominal_${1}_${2}_${3}_gammaxi${5}.root $SYSLIST $HISTLIST data_${1}_${2}_${3}_gammaxi${5}

    rm $SYSLIST
    rm $HISTLIST
done
