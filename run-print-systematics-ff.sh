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

hiBinMins=(0 20 60 100)
hiBinMaxs=(20 60 100 200)

xiBinMins=(1 0.5 1.5 3.5)
xiBinMaxs=(0 1.5 3.5 4.5)

for i1 in ${!hiBinMins[*]}
do
  hiBinMin=${hiBinMins[i1]}
  hiBinMax=${hiBinMaxs[i1]}
  for i2 in ${!xiBinMins[*]}
  do
    xiBinMin=${xiBinMins[i2]}
    xiBinMax=${xiBinMaxs[i2]}
    #echo "hiBinMin :"$hiBinMin
    #echo "hiBinMax :"$hiBinMax

    ## phoEt60, jetpt30
    touch $SYSFILELIST
    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_60_30_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_60_30_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_60_30_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_60_30_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
    ./print_systematics.exe $SYSFILELIST ff $hiBinMin $hiBinMax $xiBinMin $xiBinMax
    rm $SYSFILELIST

#    ## phoEt80, jetpt40
#    touch $SYSFILELIST
#    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_80_40_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_80_40_gxi0_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_80_40_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
#    echo -e "/export/d00/scratch/tatar/GJT-out/results/sys/data_80_40_gxi1_defnFF1-systematics.root" >> $SYSFILELIST
#    ./print_systematics.exe $SYSFILELIST ff $hiBinMin $hiBinMax $xiBinMin $xiBinMax
#    rm $SYSFILELIST
  done
done
