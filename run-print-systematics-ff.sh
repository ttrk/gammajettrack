#!/bin/bash

if [ $# -lt 0 ]; then
    echo "Usage: ./run-print-systematics-ff.sh"
    exit 1
fi

set -x

g++ print_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o print_systematics.exe || exit 1
#echo "g++ print_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o print_systematics.exe || exit 1"

SYSFILELIST=sysFileListTmp.txt
if [ -f $SYSFILELIST ]; then
        rm $SYSFILELIST
fi
echo "SYSFILELIST : $SYSFILELIST"

printRatio=0

hiBinMins=(0 20 60 100)
hiBinMaxs=(20 60 100 200)

binMins=(1 0.5 1.5 3.5)
binMaxs=(0 1.5 3.5 4.5)

doBinByBinCheck=0
if (( $doBinByBinCheck > 0 )); then

  binMins=(0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0)
  binMaxs=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5)
fi

sysDir="/home/kaya/Documents/EclipseWorkSpace/GJT/gammajettrack/results/sys/"

for i1 in ${!hiBinMins[*]}
do
  hiBinMin=${hiBinMins[i1]}
  hiBinMax=${hiBinMaxs[i1]}
  for i2 in ${!binMins[*]}
  do
    binMin=${binMins[i2]}
    binMax=${binMaxs[i2]}
    #echo "hiBinMin :"$hiBinMin
    #echo "hiBinMax :"$hiBinMax

    ## phoEt60, jetpt30
    touch $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"data_60_30_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
    ./print_systematics.exe $SYSFILELIST ff $hiBinMin $hiBinMax $binMin $binMax $printRatio
    rm $SYSFILELIST

    ## phoEt80, jetpt40
#    touch $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e $sysDir"data_80_40_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
#    ./print_systematics.exe $SYSFILELIST ff $hiBinMin $hiBinMax $binMin $binMax $printRatio
#    rm $SYSFILELIST
  done
done
