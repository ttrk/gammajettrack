#!/bin/bash

if [ $# -lt 7 ]; then
  echo "Usage: ./run-calc-systematics.sh [phoetmin] [jetptmin] [gammaxi] [sample] [type] [nominal results file] [systematics dir]"
  echo "Example: ./run-calc-systematics.sh 60 30 0 pbpbdata recoreco /export/d00/scratch/tatar/GJT-out/results/data_data_60_30_gxi0_defnFF1_ff_final.root /export/d00/scratch/tatar/GJT-out/results/sys/"
  exit 1
fi

echo "phoetmin = $1"
echo "jetptmin = $2"
echo "gammaxi  = $3"
echo "sample   = $4"
echo "type     = $5"
echo "nominal results file = $6"
echo "systematics dir = $7"

sample=$4
type=$5
nomFile=$6
varFileLR=$nomFile
sysDir=$7

g++ calc_lr_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_lr_systematics.exe || exit 1
g++ calc_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_systematics.exe || exit 1
g++ calc_ratio_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ratio_systematics.exe || exit 1
g++ calc_iso_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_iso_systematics.exe || exit 1

set -x

./calc_lr_systematics.exe $nomFile $varFileLR $4 $5 $1 $2 $3

mv longrange_${sample}_${1}_${2}_gxi${3}_defnFF1_ff_final.root $sysDir

mcsample="ppmc"
if [[ $sample == "pbpbdata" ]]; then
  mcsample="pbpbmc"
fi
./calc_iso_systematics.exe ${sysDir}nominal_iso_${mcsample}_${1}_${2}_gxi${3}_defnFF1_ff_final.root ${sysDir}iso_${mcsample}_${1}_${2}_gxi${3}_defnFF1_ff_final.root $nomFile $mcsample $sample $type $1 $2 $3
mv iso_${sample}_${1}_${2}_gxi${3}_defnFF1_ff_final.root $sysDir

## for pp : nominal = jes_qg_down
cp ${nomFile} ${sysDir}jes_qg_down_ppdata_${1}_${2}_gxi${3}_defnFF1_ff_final.root
cp ${nomFile} ${sysDir}tracking_ratio_ppdata_${1}_${2}_gxi${3}_defnFF1_ff_final.root
#####
SYSLIST=systematics_${1}_${2}_${3}_${sample}.list
if [ -f $SYSLIST ]; then
    rm $SYSLIST
fi
touch $SYSLIST

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio)

sysIndices="1 2 3 4 5 6 7 8 9 10 12 13 14"

for SYS in ${sysIndices}
do
  echo -e "${sysDir}${SYSTEMATIC[SYS]}_${sample}_${1}_${2}_gxi${3}_defnFF1_ff_final.root" >> $SYSLIST
done

HISTLIST=hist_${1}_${2}_${3}_${sample}.list
if [ -f $HISTLIST ]; then
    rm $HISTLIST
fi
touch $HISTLIST
echo -e "hff_final_${sample}_${type}_0_20" >> $HISTLIST
echo -e "hff_final_${sample}_${type}_20_60" >> $HISTLIST
echo -e "hff_final_${sample}_${type}_60_100" >> $HISTLIST
echo -e "hff_final_${sample}_${type}_100_200" >> $HISTLIST
echo -e "hff_final_${sample}_${type}_0_60" >> $HISTLIST
echo -e "hff_final_${sample}_${type}_60_200" >> $HISTLIST

cp $sysDir"/"data_${1}_${2}_gxi${3}_defnFF1-systematics.root .
./calc_systematics.exe $nomFile $SYSLIST $HISTLIST data_${1}_${2}_gxi${3}_defnFF1
mv sys_hff_final_${sample}_${type}_*_*-data_${1}_${2}_gxi${3}_defnFF1.png $sysDir
mv data_${1}_${2}_gxi${3}_defnFF1-systematics.root $sysDir

rm $SYSLIST
rm $HISTLIST

echo "running ratio systematics"
SYSHISTLIST=syshist_${1}_${2}_${3}.list
if [ -f $SYSHISTLIST ]; then
  rm $SYSHISTLIST
fi
echo -e "0_20" >> $SYSHISTLIST
echo -e "20_60" >> $SYSHISTLIST
echo -e "60_100" >> $SYSHISTLIST
echo -e "100_200" >> $SYSHISTLIST
echo -e "0_60" >> $SYSHISTLIST
echo -e "60_200" >> $SYSHISTLIST

SYSFILELIST=sysfile_${1}_${2}_${3}.list
if [ -f $SYSFILELIST ]; then
  rm $SYSFILELIST
fi
echo -e "data_${1}_${2}_gxi${3}_defnFF1-systematics.root" >> $SYSFILELIST
echo -e "data_${1}_${2}_gxi${3}_defnFF1-systematics.root" >> $SYSFILELIST

cp $sysDir"/"data_${1}_${2}_gxi${3}_defnFF1-systematics.root .
./calc_ratio_systematics.exe ff $SYSFILELIST $SYSHISTLIST data_${1}_${2}_gxi${3}_defnFF1
mv sys_hff_final_*_*_*-data_${1}_${2}_gxi${3}_defnFF1.png $sysDir
mv data_${1}_${2}_gxi${3}_defnFF1-systematics.root $sysDir

rm $SYSHISTLIST
rm $SYSFILELIST





