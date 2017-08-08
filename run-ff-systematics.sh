#!/bin/bash

if [ $# -lt 7 ]; then
  echo "Usage: ./run-ff-systematics.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [nominal results]"
  echo "Example: ./run-ff-systematics.sh 60 1000 30 1 0 pbpbdata [nominal results]"
  exit 1
fi

echo "phoetmin = $1"
echo "phoetmax = $2"
echo "jetptmin = $3"
echo "trkptmin = $4"
echo "gammaxi  = $5"
echo "sample   = $6"
echo "nominal results = $7"

kFF=1 # 0 : old definition, 1 : new definition

if [ $6 = "pbpbdata" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-Data-skim-170712.root"
    TYPE="recoreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170712.root"
    MCSAMPLE="pbpbmc"
elif [ $6 = "ppdata" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/pp-Data-skim-170712.root"
    TYPE="srecoreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/pp-MC-skim-170712.root"
    MCSAMPLE="ppmc"
else
    echo "invalid sample"
    exit 1
fi

nomDir="/export/d00/scratch/"$USER"/GJT-out/results/"
sysDir="/export/d00/scratch/"$USER"/GJT-out/results/sys/"
mkdir -p $nomDir
mkdir -p $sysDir
echo "directory for nominal results : $nomDir"
echo "directory for sys variations  : $sysDir"
nomPrefix=$nomDir"data"
sysPrefix=$sysDir

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down longrange)

echo "compiling macros..."
g++ jetff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetff.exe || exit 1
g++ draw_ff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff.exe || exit 1
g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1
g++ calc_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_systematics.exe || exit 1
g++ calc_ratio_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ratio_systematics || exit 1
g++ calc_iso_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_iso_systematics.exe || exit 1

set -x

for SYS in 1 2 3 4 6 9 10
do
    ./jetff.exe $SKIM $6 0 20 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
    ./jetff.exe $SKIM $6 20 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
    ./jetff.exe $SKIM $6 60 100 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
    ./jetff.exe $SKIM $6 100 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
    ./jetff.exe $SKIM $6 0 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
    ./jetff.exe $SKIM $6 60 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
done
wait

./draw_ff.exe $6 ${nomPrefix}_data_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_up_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 2 $TYPE
./draw_ff.exe $6 ${nomPrefix}_data_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_down_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} -2 $TYPE

# isolation systematics
./jetff.exe $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 0 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 60 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}nominal_iso 0 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 0 20 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 20 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 60 100 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 100 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 0 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
./jetff.exe $MCSKIM $MCSAMPLE 60 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}iso 5 $kFF &
wait

hadd -f ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
rm ${sysPrefix}nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
hadd -f ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
./draw_ff.exe $MCSAMPLE ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}

hadd -f ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
rm ${sysPrefix}iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
hadd -f ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
./draw_ff.exe $MCSAMPLE ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}

./calc_iso_systematics.exe ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root $7 $MCSAMPLE $6 $TYPE $1 $3 $5
mv iso_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root $sysDir

for SYS in 1 2 3 4 6 9 10
do
    hadd -f ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
    rm ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
    hadd -f ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
    ./draw_ff.exe $6 ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}
done

# SYSLIST=systematics_${1}_${3}_${5}_${6}.list
# if [ -f $SYSLIST ]; then
#     rm $SYSLIST
# fi
# touch $SYSLIST

# for SYS in 1 2 3 4 5 6 7 8 9 10
# do
#     echo -e "${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root" >> $SYSLIST
# done

# # placeholder for systematics yet to be implemented
# # none

# HISTLIST=hist_${1}_${3}_${5}_${6}.list
# if [ -f $HISTLIST ]; then
#     rm $HISTLIST
# fi
# touch $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_0_20" >> $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_20_60" >> $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_60_100" >> $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_100_200" >> $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_0_60" >> $HISTLIST
# echo -e "hff_final_${6}_${TYPE}_60_200" >> $HISTLIST

# cp $sysDir"/"data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root .
# ./calc_systematics.exe $7 $SYSLIST $HISTLIST data_${1}_${3}_gxi${5}_defnFF${kFF}
# mv sys_hff_final_${6}_${TYPE}_*_*-data_${1}_${3}_gxi${5}_defnFF${kFF}.png $sysDir
# mv data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root $sysDir

# rm $SYSLIST
# rm $HISTLIST

# echo "running ratio systematics"
# SYSHISTLIST=syshist_${1}_${3}_${5}.list
# if [ -f $SYSHISTLIST ]; then
#   rm $SYSHISTLIST
# fi
# echo -e "0_20" >> $SYSHISTLIST
# echo -e "20_60" >> $SYSHISTLIST
# echo -e "60_100" >> $SYSHISTLIST
# echo -e "100_200" >> $SYSHISTLIST
# echo -e "0_60" >> $SYSHISTLIST
# echo -e "60_200" >> $SYSHISTLIST

# SYSFILELIST=sysfile_${1}_${3}_${5}.list
# if [ -f $SYSFILELIST ]; then
#   rm $SYSFILELIST
# fi
# echo -e "data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root" >> $SYSFILELIST
# echo -e "data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root" >> $SYSFILELIST

# cp $sysDir"/"data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root .
# ./calc_ratio_systematics ff $SYSFILELIST $SYSHISTLIST data_${1}_${3}_gxi${5}_defnFF${kFF}
# mv sys_hff_final_*_*_*-data_${1}_${3}_gxi${5}_defnFF${kFF}.png $sysDir
# mv data_${1}_${3}_gxi${5}_defnFF${kFF}-systematics.root $sysDir

# rm $SYSHISTLIST
# rm $SYSFILELIST



