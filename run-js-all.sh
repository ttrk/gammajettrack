#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "Usage: ./run-js-all.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi]"
    echo "Example: ./run-js-all.sh 60 1000 30 1 0"
    exit 1
fi

echo "compiling macros..."
make jetshape calc_systematics calc_ratio_systematics plot_js

set -x

echo "running data"
./run-js-data.sh $@ both nominal

echo "running systematics"
./run-js-systematics.sh $@ pbpbdata nominal_pbpbdata_${1}_${3}_gxi${5}_js_final.root
./run-js-systematics.sh $@ ppdata nominal_ppdata_${1}_${3}_gxi${5}_js_final.root

echo "running ratio systematics"
SYSHISTLIST=syshist_${1}_${3}_${5}.list
[ -f $SYSHISTLIST ] && rm $SYSHISTLIST

touch $SYSHISTLIST
echo -e "0_20" >> $SYSHISTLIST
echo -e "20_60" >> $SYSHISTLIST
echo -e "60_100" >> $SYSHISTLIST
echo -e "100_200" >> $SYSHISTLIST

SYSFILELIST=sysfile_${1}_${3}_${5}.list
if [ -f $SYSFILELIST ]; then
    rm $SYSFILELIST
fi
touch $SYSFILELIST

echo -e "data_${1}_${3}_gxi${5}-systematics.root" >> $SYSFILELIST
echo -e "data_${1}_${3}_gxi${5}-systematics.root" >> $SYSFILELIST

./calc_ratio_systematics js $SYSFILELIST $SYSHISTLIST data_${1}_${3}_gxi${5}

rm $SYSHISTLIST
rm $SYSFILELIST

DATAFILE=nominal_data_${1}_${3}_gxi${5}_js_final.root
SYSFILE=data_${1}_${3}_gxi${5}-systematics.root

echo "plotting final results"
PLOTLIST=plot_${1}_${3}_${5}_final.list
[ -f $PLOTLIST ] && rm $PLOTLIST

touch $PLOTLIST
echo -e "pp (smeared)" >> $PLOTLIST
echo -e "hjetshape_final_ppdata_srecoreco_0_20" >> $PLOTLIST
echo -e "hjetshape_final_ppdata_srecoreco_20_60" >> $PLOTLIST
echo -e "hjetshape_final_ppdata_srecoreco_60_100" >> $PLOTLIST
echo -e "hjetshape_final_ppdata_srecoreco_100_200" >> $PLOTLIST
echo -e "PbPb" >> $PLOTLIST
echo -e "hjetshape_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
echo -e "hjetshape_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
echo -e "hjetshape_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
echo -e "hjetshape_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

./plot_js nominal_data_${1}_${3}_gxi${5}_js_final.root final_js_${1}_${3}_gxi${5} $PLOTLIST 1 $5 $1 $3 0 data_${1}_${3}_gxi${5}-systematics.root

rm $PLOTLIST
