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
    SKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
    TYPE="recoreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
    MCSAMPLE="pbpbmc"
elif [ $6 = "pbpbmc" ]; then
    SKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
    TYPE="recoreco"
    MCSKIM="/export/d00/scratch/biran/photon-jet-track/PbPb-MC-skim-170911.root"
    MCSAMPLE="pbpbmc"
elif [ $6 = "ppdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/pp-Data-skim-170911.root"
    TYPE="srecoreco"
    MCSKIM="/export/d00/scratch/tatar/GJT-out/pp-MC-skim-170911.root"
    MCSAMPLE="ppmc"
elif [ $6 = "ppmc" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/pp-MC-skim-170911.root"
    TYPE="srecoreco"
    MCSKIM="/export/d00/scratch/tatar/GJT-out/pp-MC-skim-170911.root"
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

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio phoeffcorr)

echo "compiling macros..."
g++ jetff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetff.exe || exit 1
g++ draw_ff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff.exe || exit 1
g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1

set -x

sysIndices="1 2 3 4 6 9 10 15"
runSysIso=1
if [ $6 = "pbpbdata" ] || [ $6 = "pbpbmc" ]; then
  sysIndices="1 2 3 4 6 9 10 11 12 14 15"
  runSysIso=1
fi

for SYS in ${sysIndices}
do
  ./jetff.exe $SKIM $6 0 20 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
  ./jetff.exe $SKIM $6 20 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
  ./jetff.exe $SKIM $6 60 100 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
  ./jetff.exe $SKIM $6 100 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
  ./jetff.exe $SKIM $6 0 60 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
  ./jetff.exe $SKIM $6 60 200 $1 $2 $3 $TYPE $4 $5 ${sysPrefix}${SYSTEMATIC[SYS]} $SYS $kFF &
done
wait

if [ $runSysIso == 1 ]; then
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
fi

if [ $6 = "pbpbdata" ] || [ $6 = "ppdata" ]; then
  ./draw_ff.exe $6 ${nomPrefix}_data_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_up_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 2 $TYPE
  ./draw_ff.exe $6 ${nomPrefix}_data_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_down_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} -2 $TYPE
elif [ $6 = "pbpbmc" ] || [ $6 = "ppmc" ]; then
  ./draw_ff.exe $6 ${nomDir}"closure/ffclosure_"${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_up_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 2 $TYPE
  ./draw_ff.exe $6 ${nomDir}"closure/ffclosure_"${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}purity_down_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} -2 $TYPE
fi

for SYS in ${sysIndices}
do
  hadd -f ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  rm ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  hadd -f ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
  ./draw_ff.exe $6 ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}${SYSTEMATIC[SYS]}_${6}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}
done

if [ $runSysIso == 1 ]; then
  hadd -f ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  rm ${sysPrefix}nominal_iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  hadd -f ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
  ./draw_ff.exe $MCSAMPLE ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}nominal_iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}

  hadd -f ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root ${sysPrefix}iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  rm ${sysPrefix}iso_${MCSAMPLE}_${TYPE}_${1}_${3}_${5}_${kFF}_*_*.root
  hadd -f ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_${TYPE}_ff.root
  ./draw_ff.exe $MCSAMPLE ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_merged.root ${sysPrefix}iso_${MCSAMPLE}_${1}_${3}_gxi${5}_defnFF${kFF}_ff_final.root ${1} 0 ${TYPE}
fi

