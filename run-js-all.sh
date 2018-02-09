#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "Usage: ./run-js-all.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi]"
    echo "Example: ./run-js-all.sh 60 1000 30 1 0"
    exit 1
fi

echo "compiling macros..."
make calc_js_ratio_systematics plot_js

DATAFILE=nominal_data_${1}_${3}_gxi${5}_js_final.root
SYSFILE=nominal_data_${1}_${3}_gxi${5}-systematics.root

set -x

echo "running data"
./run-js-data.sh $@ both nominal

echo "running systematics"
./run-js-systematics.sh $@ pbpbdata nominal_pbpbdata_${1}_${3}_gxi${5}_js_final.root
./run-js-systematics.sh $@ ppdata nominal_ppdata_${1}_${3}_gxi${5}_js_final.root

hadd -f $SYSFILE nominal_ppdata_${1}_${3}_gxi${5}-systematics.root nominal_pbpbdata_${1}_${3}_gxi${5}-systematics.root

echo "running ratio systematics"
SYSHISTLIST=syshist_${1}_${3}_${5}.list
[ -f $SYSHISTLIST ] && rm $SYSHISTLIST

touch $SYSHISTLIST
echo -e "0_20" >> $SYSHISTLIST
echo -e "20_60" >> $SYSHISTLIST
echo -e "60_100" >> $SYSHISTLIST
echo -e "100_200" >> $SYSHISTLIST

PBPBLIST=systematics_${1}_${3}_${5}_pbpbdata.list
PPLIST=systematics_${1}_${3}_${5}_ppdata.list

./calc_js_ratio_systematics $DATAFILE $PBPBLIST $PPLIST $SYSHISTLIST nominal_data_${1}_${3}_gxi${5}

rm $SYSHISTLIST $PBPBLIST $PPLIST

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

./plot_js $DATAFILE final_js_${1}_${3}_gxi${5} $PLOTLIST 1 $5 $1 $3 0 $SYSFILE

rm $PLOTLIST
