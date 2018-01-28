#!/bin/bash

if [ $# -lt 9 ]; then
    echo "Usage: ./quick-js-plot.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [hist label] [types...]"
    echo "Example: ./quick-js-plot.sh 60 1000 30 1 0 pbpbmc closure final sgengen sgenreco recogen recoreco"
    exit 1
fi

echo "compiling macros..."
make plot_js

arglist=""
for i in ${@:9}; do
    arglist+="_"
    arglist+="$i"
done

PLOTLIST=jsplot_${1}_${3}_${5}_${7}_${8}${arglist}.list
if [ -f $PLOTLIST ]; then
    rm $PLOTLIST
fi
touch $PLOTLIST

set -x

if [ $9 = "data" ]; then
    echo -e "pp (smeared)" >> $PLOTLIST
    echo -e "hjetshape_${8}_ppdata_srecoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_${8}_ppdata_srecoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_${8}_ppdata_srecoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_${8}_ppdata_srecoreco_100_200" >> $PLOTLIST
    echo -e "PbPb" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_js ${7}_data_${1}_${3}_gxi${5}_js_final.root ${7}_data_gxi${5}_${1}_${3}_${8} $PLOTLIST 1 $5 $1 $3
elif [ $9 = "datamc" ]; then
    echo -e "pp MC" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppmc_recoreco_0_20" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppmc_recoreco_20_60" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppmc_recoreco_60_100" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppmc_recoreco_100_200" >> ${PLOTLIST}
    echo -e "pp data" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppdata_recoreco_0_20" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppdata_recoreco_20_60" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppdata_recoreco_60_100" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_ppdata_recoreco_100_200" >> ${PLOTLIST}

    ./plot_js ${7}_datamc_${1}_${3}_gxi${5}_js_final.root ${7}_datamc_gxi${5}_${1}_${3}_${8} ${PLOTLIST} 1 $5 $1 $3
    rm $PLOTLIST

    echo -e "PbPb MC" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbmc_recoreco_0_20" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbmc_recoreco_20_60" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbmc_recoreco_60_100" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbmc_recoreco_100_200" >> ${PLOTLIST}
    echo -e "PbPb data" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbdata_recoreco_0_20" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbdata_recoreco_20_60" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbdata_recoreco_60_100" >> ${PLOTLIST}
    echo -e "hjetshape_${8}_pbpbdata_recoreco_100_200" >> ${PLOTLIST}

    ./plot_js ${7}_datamc_${1}_${3}_gxi${5}_js_final.root ${7}_datamc_gxi${5}_${1}_${3}_${8} ${PLOTLIST} 1 $5 $1 $3
elif [ $9 = "pbpbpp" ]; then
    hadd -f ${7}_pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pbpb${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pp${6}_${1}_${3}_gxi${5}_js_final.root

    echo -e "pp ${11}" >> $PLOTLIST
    echo -e "hjetshape_${8}_pp${6}_${11}_0_20" >> $PLOTLIST
    echo -e "hjetshape_${8}_pp${6}_${11}_20_60" >> $PLOTLIST
    echo -e "hjetshape_${8}_pp${6}_${11}_60_100" >> $PLOTLIST
    echo -e "hjetshape_${8}_pp${6}_${11}_100_200" >> $PLOTLIST
    echo -e "pbpb ${10}" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpb${6}_${10}_0_20" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpb${6}_${10}_20_60" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpb${6}_${10}_60_100" >> $PLOTLIST
    echo -e "hjetshape_${8}_pbpb${6}_${10}_100_200" >> $PLOTLIST

    ./plot_js ${7}_pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_pbpb_${10}_pp_${11}_gxi${5}_${1}_${3}_${8} $PLOTLIST 1 $5 $1 $3
else
    for i in ${@:9}; do
        echo -e "$i" >> $PLOTLIST
        echo -e "hjetshape_${8}_${6}_${i}_0_20" >> $PLOTLIST
        echo -e "hjetshape_${8}_${6}_${i}_20_60" >> $PLOTLIST
        echo -e "hjetshape_${8}_${6}_${i}_60_100" >> $PLOTLIST
        echo -e "hjetshape_${8}_${6}_${i}_100_200" >> $PLOTLIST
    done

    [[ $6 = "pp"* ]] && opt=1

    ./plot_js ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_gxi${5}_${1}_${3}_${8}${arglist} $PLOTLIST 1 $5 $1 $3 $opt
fi

rm $PLOTLIST
