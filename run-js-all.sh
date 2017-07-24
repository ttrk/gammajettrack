#!/bin/bash

if [[ $# -lt 5 ]]; then
    echo "Usage: ./run-js-all.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [[final plots only]]"
    echo "Example: ./run-js-all.sh 80 1000 40 1 0"
    exit 1
fi

echo "compiling macros..."
# g++ jetshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetshape || exit 1
# g++ draw_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_js || exit 1
# g++ calc_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_systematics || exit 1
# g++ calc_ratio_systematics.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ratio_systematics || exit 1
# g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results || exit 1

set -x

if [[ $6 -ne 1 ]]; then
    echo "running data"
    ./run-js-data.sh $@ 0

    echo "running systematics"
    ./run-js-systematics.sh $@ pbpbdata data_pbpbdata_${1}_${3}_gxi${5}_js_final.root
    ./run-js-systematics.sh $@ ppdata data_ppdata_${1}_${3}_gxi${5}_js_final.root

    echo "running ratio systematics"
    SYSHISTLIST=syshist_${1}_${3}_${5}.list
    if [ -f $SYSHISTLIST ]; then
        rm $SYSHISTLIST
    fi
    echo -e "0_20" >> $SYSHISTLIST
    echo -e "20_60" >> $SYSHISTLIST
    echo -e "60_100" >> $SYSHISTLIST
    echo -e "100_200" >> $SYSHISTLIST

    SYSFILELIST=sysfile_${1}_${3}_${5}.list
    if [ -f $SYSFILELIST ]; then
        rm $SYSFILELIST
    fi
    echo -e "data_${1}_${3}_gxi${5}-systematics.root" >> $SYSFILELIST
    echo -e "data_${1}_${3}_gxi${5}-systematics.root" >> $SYSFILELIST

    ./calc_ratio_systematics js $SYSFILELIST $SYSHISTLIST data_${1}_${3}_gxi${5}

    rm $SYSHISTLIST
    rm $SYSFILELIST
fi

DATAFILE=data_data_${1}_${3}_gxi${5}_js_final.root
SYSFILE=data_${1}_${3}_gxi${5}-systematics.root

if [[ $6 -eq 1 ]]; then
    if [ ! -f $DATAFILE ]; then
        echo "data file: $DATAFILE not found!"
        exit 1
    fi
    if [ ! -f $SYSFILE ]; then
        echo "systematics file: $SYSFILE not found!"
        exit 1
    fi
fi

echo "plotting final results"
PLOTLIST=plot_${1}_${3}_${5}_final.list
if [ -f $PLOTLIST ]; then
    rm $PLOTLIST
fi
touch $PLOTLIST

echo -e "pp (smeared)" >> $PLOTLIST
echo -e "hjs_final_ppdata_srecoreco_0_20" >> $PLOTLIST
echo -e "hjs_final_ppdata_srecoreco_20_60" >> $PLOTLIST
echo -e "hjs_final_ppdata_srecoreco_60_100" >> $PLOTLIST
echo -e "hjs_final_ppdata_srecoreco_100_200" >> $PLOTLIST
echo -e "PbPb" >> $PLOTLIST
echo -e "hjs_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
echo -e "hjs_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
echo -e "hjs_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
echo -e "hjs_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

./plot_results data_data_${1}_${3}_gxi${5}_js_final.root final_js_${1}_${3}_gxi${5} $PLOTLIST 1 $5 $1 $3 0 data_${1}_${3}_gxi${5}-systematics.root

rm $PLOTLIST
