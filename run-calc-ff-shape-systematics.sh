#!/bin/bash

if [ $# -lt 7 ]; then
  echo "Usage: ./run-calc-ff-shape-systematics.sh [phoetmin] [jetptmin] [gammaxi] [obs] [sample] [type] [nominal results file] [systematics dir]"
  echo "Example: ./run-calc-ff-shape-systematics.sh 60 30 0 2 pbpbdata recoreco /export/d00/scratch/tatar/GJT-out/results/jsdata_pbpbdata_60_30_gxi0_obs2_ffjs_final.root /export/d00/scratch/tatar/GJT-out/results/sys/"
  exit 1
fi

phoetmin=$1
jetptmin=$2
gammaxi=$3
obs=$4
sample=$5
type=$6
nomFile=$7
sysDir=$8

echo "phoetmin = $phoetmin"
echo "jetptmin = $jetptmin"
echo "gammaxi  = $gammaxi"
echo "obs      = $obs"
echo "sample   = $sample"
echo "type     = $type"
echo "nominal results file = $nomFile"
echo "systematics dir = ${sysDir}"

g++ calc_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_systematics.exe || exit 1
g++ calc_ratio_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ratio_systematics.exe || exit 1
g++ calc_iso_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_iso_systematics.exe || exit 1


obsSTR="js"
label="jssys"
histPrefix="hjs"
if [ $obs = "1" ]; then
  obsSTR="ff"
  label="ffsys"
  histPrefix="hff"
fi

set -x

if [ $obs = "1" ]; then
  g++ calc_lr_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_lr_systematics.exe || exit 1
  varFileLR=$sysDir/${label}_longrange_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  outputTMP=longrange_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_defnFF1_ff_final.root
  cp $varFileLR $outputTMP
  ./calc_lr_systematics.exe $nomFile $varFileLR $sample $type $phoetmin $jetptmin $gammaxi hff
  mv $outputTMP $varFileLR
fi

mcsample="ppmc"
if [ $sample = "pbpbdata" ]; then
  mcsample="pbpbmc"
fi

inputNominalIso=${sysDir}/${label}_nominal_iso_${mcsample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
inputIso=${sysDir}/${label}_iso_${mcsample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
./calc_iso_systematics.exe $inputNominalIso $inputIso $nomFile $mcsample $sample $type $phoetmin $jetptmin $gammaxi ${label}_iso
outputTmp=${label}_iso_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_js_final.root
if [ $obs = "1" ]; then
  outputTmp=iso_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_defnFF1_ff_final.root
fi
outputIso=${sysDir}/${label}_iso_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
mv $outputTmp $outputIso

## for pp : nominal = jes_qg_down
if [ $sample = "ppdata" ]; then
  cp ${nomFile} ${sysDir}/${label}_jes_qg_down_ppdata_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  cp ${nomFile} ${sysDir}/${label}_tracking_ratio_ppdata_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  cp ${nomFile} ${sysDir}/${label}_noUEscale_ppdata_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
fi
#####
SYSLIST=systematics_${phoetmin}_${jetptmin}_${gammaxi}_${sample}.list
if [ -f $SYSLIST ]; then
    rm $SYSLIST
fi
touch $SYSLIST

#SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio eta_reflection etagt0p3)
SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio phoeffcorr noUEscale)

sysIndices="1 2 3 4 5 6 7 8 9 10 12 13 14 15 16"

for SYS in ${sysIndices}
do
  echo -e "${sysDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root" >> $SYSLIST
done

HISTLIST=hist_${phoetmin}_${jetptmin}_${gammaxi}_${sample}.list
if [ -f $HISTLIST ]; then
    rm $HISTLIST
fi
touch $HISTLIST
echo -e "${histPrefix}_final_${sample}_${type}_100_200" >> $HISTLIST
if [ $sample = "pbpbdata" ] || ([ $sample = "ppdata" ] && ([ $type = "srecoreco" ] || [ $type = "scorrjsrecoreco" ])); then
  echo -e "${histPrefix}_final_${sample}_${type}_60_100" >> $HISTLIST
  echo -e "${histPrefix}_final_${sample}_${type}_20_60" >> $HISTLIST
  echo -e "${histPrefix}_final_${sample}_${type}_0_20" >> $HISTLIST
  echo -e "${histPrefix}_final_${sample}_${type}_60_200" >> $HISTLIST
  echo -e "${histPrefix}_final_${sample}_${type}_0_60" >> $HISTLIST
fi

cp ${sysDir}"/"${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root .
./calc_systematics.exe $nomFile $SYSLIST $HISTLIST ${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs
mv sys_${histPrefix}_final_${sample}_${type}_*_*-data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs.png ${sysDir}
mv ${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root ${sysDir}

rm $SYSLIST
rm $HISTLIST

echo "running ratio systematics"
SYSHISTLIST=syshist_${phoetmin}_${jetptmin}_${gammaxi}.list
if [ -f $SYSHISTLIST ]; then
  rm $SYSHISTLIST
fi
echo -e "100_200" >> $SYSHISTLIST
echo -e "60_100" >> $SYSHISTLIST
echo -e "20_60" >> $SYSHISTLIST
echo -e "0_20" >> $SYSHISTLIST
echo -e "60_200" >> $SYSHISTLIST
echo -e "0_60" >> $SYSHISTLIST

SYSFILELIST=sysfile_${phoetmin}_${jetptmin}_${gammaxi}.list
if [ -f $SYSFILELIST ]; then
  rm $SYSFILELIST
fi
echo -e "${sysDir}/${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root" >> $SYSFILELIST
echo -e "${sysDir}/${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root" >> $SYSFILELIST

cp ${sysDir}"/"${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root .
./calc_ratio_systematics.exe $obsSTR $SYSFILELIST $SYSHISTLIST ${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs
mv sys_${histPrefix}_final_*_*_*-${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs.png ${sysDir}
mv ${label}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root ${sysDir}

rm $SYSHISTLIST
rm $SYSFILELIST


