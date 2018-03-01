#!/bin/bash

if [ $# -lt 0 ]; then
    echo "Usage: ./run-print-systematics-js.sh"
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

binMins=(1 0.0 0.1 0.2)
binMaxs=(0 0.1 0.2 0.3)
#binMins=(0.00 0.05 0.10 0.15 0.20 0.25)
#binMaxs=(0.05 0.10 0.15 0.20 0.25 0.30)


sysDir="/export/d00/scratch/tatar/GJT-out/results/sys/"

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
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    echo -e $sysDir"jssys_data_60_30_gxi0_obs2_ffjs-systematics.root" >> $SYSFILELIST
    ./print_systematics.exe $SYSFILELIST js $hiBinMin $hiBinMax $binMin $binMax 1
    rm $SYSFILELIST

  done
done
