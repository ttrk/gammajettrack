#!/bin/bash

if [ $# -lt 7 ]; then
  echo "Usage: ./run-js-systematics.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [nominal results]"
  echo "Example: ./run-js-systematics.sh 60 1000 30 1 0 pbpbdata [nominal results]"
  exit 1
fi

if [ $6 = "pbpbdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
    TYPE="recoreco"
    GENTYPE="gengen"
    RECOTYPE="genreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
    MCSAMPLE="pbpbmc"
elif [ $6 = "ppdata" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170911.root"
    TYPE="srecoreco"
    GENTYPE="gengen"
    RECOTYPE="genreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-170911.root"
    MCSAMPLE="ppmc"
else
    echo "invalid sample"
    exit 1
fi

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking placeholder jes_gluon jes_quark longrange etareflect)

echo "compiling macros..."
make jetshape draw_js calc_iso_systematics calc_js_systematics

set -x

for SYS in 1 2 3 4 6 11 12 13 14
do
    ./jetshape $SKIM $6 0 20 $1 $2 $3 $TYPE $4 $5 ${SYSTEMATIC[SYS]} $SYS &
    ./jetshape $SKIM $6 20 60 $1 $2 $3 $TYPE $4 $5 ${SYSTEMATIC[SYS]} $SYS &
    ./jetshape $SKIM $6 60 100 $1 $2 $3 $TYPE $4 $5 ${SYSTEMATIC[SYS]} $SYS &
    ./jetshape $SKIM $6 100 200 $1 $2 $3 $TYPE $4 $5 ${SYSTEMATIC[SYS]} $SYS &
done
wait

./draw_js $6 nominal_data_${1}_${3}_gxi${5}_js_merged.root purity_up_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 2 $TYPE
./draw_js $6 nominal_data_${1}_${3}_gxi${5}_js_merged.root purity_down_${6}_${1}_${3}_gxi${5}_js_final.root ${1} -2 $TYPE

# isolation systematics
./jetshape $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $TYPE $4 $5 nominal_iso 0 &
./jetshape $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $TYPE $4 $5 nominal_iso 0 &
./jetshape $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $TYPE $4 $5 nominal_iso 0 &
./jetshape $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $TYPE $4 $5 nominal_iso 0 &
./jetshape $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $TYPE $4 $5 iso 5 &
./jetshape $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $TYPE $4 $5 iso 5 &
./jetshape $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $TYPE $4 $5 iso 5 &
./jetshape $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $TYPE $4 $5 iso 5 &
wait

hadd -f nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_${TYPE}_js.root nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_*_*.root
rm nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_*_*.root
hadd -f nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_${TYPE}_js.root
./draw_js $MCSAMPLE nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${TYPE}

hadd -f iso_${MCSAMPLE}_${1}_${3}_gxi${5}_${TYPE}_js.root iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_*_*.root
rm iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_*_*.root
hadd -f iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root iso_${MCSAMPLE}_${1}_${3}_gxi${5}_${TYPE}_js.root
./draw_js $MCSAMPLE iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${TYPE}

./calc_iso_systematics nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root iso_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root $7 $MCSAMPLE $6 $TYPE $1 $3 $5

# tracking systematics
./jetshape $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $GENTYPE $4 $5 trk_gen &
./jetshape $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $GENTYPE $4 $5 trk_gen &
./jetshape $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $GENTYPE $4 $5 trk_gen &
./jetshape $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $GENTYPE $4 $5 trk_gen &
./jetshape $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $RECOTYPE $4 $5 trk_reco &
./jetshape $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $RECOTYPE $4 $5 trk_reco &
./jetshape $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $RECOTYPE $4 $5 trk_reco &
./jetshape $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $RECOTYPE $4 $5 trk_reco &
wait

hadd -f trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_${GENTYPE}_js.root trk_gen_${MCSAMPLE}_${GENTYPE}_${1}_${3}_${5}_*_*.root
rm trk_gen_${MCSAMPLE}_${GENTYPE}_${1}_${3}_${5}_*_*.root
hadd -f trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_${GENTYPE}_js.root
./draw_js $MCSAMPLE trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${GENTYPE}

hadd -f trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_${RECOTYPE}_js.root trk_reco_${MCSAMPLE}_${RECOTYPE}_${1}_${3}_${5}_*_*.root
rm trk_reco_${MCSAMPLE}_${RECOTYPE}_${1}_${3}_${5}_*_*.root
hadd -f trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_${RECOTYPE}_js.root
./draw_js $MCSAMPLE trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_js_merged.root trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${RECOTYPE}

./calc_iso_systematics trk_gen_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root trk_reco_${MCSAMPLE}_${1}_${3}_gxi${5}_js_final.root $7 $MCSAMPLE $6 $TYPE $1 $3 $5 tracking $GENTYPE $RECOTYPE

for SYS in 1 2 3 4 6 11 12 13 14
do
    hadd -f ${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_${TYPE}_js.root ${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_*_*.root
    rm ${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_*_*.root
    hadd -f ${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_js_merged.root ${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_${TYPE}_js.root
    ./draw_js $6 ${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_js_merged.root ${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${TYPE}
done

SYSLIST=systematics_${1}_${3}_${5}_${6}.list
[ -f $SYSLIST ] && rm $SYSLIST

touch $SYSLIST
for SYS in 1 2 3 4 5 6 7 8 9 11 12 13 14
do
    echo -e "${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_js_final.root" >> $SYSLIST
done

# placeholder for systematics yet to be implemented
# none

HISTLIST=hist_${1}_${3}_${5}_${6}.list
[ -f $HISTLIST ] && rm $HISTLIST

touch $HISTLIST
echo -e "hjetshape_final_${6}_${TYPE}_0_20" >> $HISTLIST
echo -e "hjetshape_final_${6}_${TYPE}_20_60" >> $HISTLIST
echo -e "hjetshape_final_${6}_${TYPE}_60_100" >> $HISTLIST
echo -e "hjetshape_final_${6}_${TYPE}_100_200" >> $HISTLIST

./calc_js_systematics $7 $SYSLIST $HISTLIST nominal_${6}_${1}_${3}_gxi${5}

rm $HISTLIST
