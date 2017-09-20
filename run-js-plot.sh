#!/bin/bash

if [ $# -lt 8 ]; then
  echo "Usage: ./run-js-plot.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
  echo "Example: ./run-js-plot.sh 60 1000 30 1 0 pbpbmc closure sgengen sgenreco recogen recoreco"
  exit 1
fi

echo "compiling macros..."
make draw_js plot_results

PLOTLIST=jsplot_${1}_${3}_${5}_${6}_${7}_${8}.list
if [ -f $PLOTLIST ]; then
    rm $PLOTLIST
fi
touch $PLOTLIST

set -x

if [ $8 = "data" ]; then
    hadd -f data_data_${1}_${3}_gxi${5}_js_merged.root data_ppdata_${1}_${3}_gxi${5}_srecoreco_js.root data_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root
    ./draw_js pbpbdata data_data_${1}_${3}_gxi${5}_js_merged.root data_data_${1}_${3}_gxi${5}_js_final.root ${1} 0 recoreco
    ./draw_js ppdata data_data_${1}_${3}_gxi${5}_js_merged.root data_data_${1}_${3}_gxi${5}_js_final.root ${1} 0 srecoreco

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

    ./plot_results data_data_${1}_${3}_gxi${5}_js_final.root data_data_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
elif [ $8 = "pbpbpp" ]; then
    hadd -f ${7}_pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pbpb${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pp${6}_${1}_${3}_gxi${5}_js_final.root
    echo -e "pp ${10}" >> $PLOTLIST
    echo -e "hjetshape_final_pp${6}_${10}_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pp${6}_${10}_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pp${6}_${10}_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pp${6}_${10}_100_200" >> $PLOTLIST
    echo -e "pbpb ${9}" >> $PLOTLIST
    echo -e "hjetshape_final_pbpb${6}_${9}_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pbpb${6}_${9}_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pbpb${6}_${9}_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pbpb${6}_${9}_100_200" >> $PLOTLIST

    ./plot_results ${7}_pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_pbpb_${9}_pp_${10}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
elif [ $8 = "datamc" ]; then
    hadd -f datamc_${1}_${3}_gxi${5}_js_merged.root data_ppdata_${1}_${3}_gxi${5}_recoreco_js.root data_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root closure_ppmc_${1}_${3}_gxi${5}_recoreco_js.root closure_pbpbmc_${1}_${3}_gxi${5}_recoreco_js.root
    ./draw_js pbpbdata datamc_${1}_${3}_gxi${5}_js_merged.root datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppdata datamc_${1}_${3}_gxi${5}_js_merged.root datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js pbpbmc datamc_${1}_${3}_gxi${5}_js_merged.root datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppmc datamc_${1}_${3}_gxi${5}_js_merged.root datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco

    echo -e "pp MC" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppmc_recoreco_0_20" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppmc_recoreco_20_60" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppmc_recoreco_60_100" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppmc_recoreco_100_200" >> ${PLOTLIST}_pp
    echo -e "pp data" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppdata_recoreco_0_20" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppdata_recoreco_20_60" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppdata_recoreco_60_100" >> ${PLOTLIST}_pp
    echo -e "hjetshape_final_ppdata_recoreco_100_200" >> ${PLOTLIST}_pp

    ./plot_results datamc_${1}_${3}_gxi${5}_js_final.root datamc_gxi${5}_${1}_${3}_pp ${PLOTLIST}_pp 1 $5 $1 $3

    echo -e "PbPb MC" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbmc_recoreco_0_20" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbmc_recoreco_20_60" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbmc_recoreco_60_100" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbmc_recoreco_100_200" >> ${PLOTLIST}_pbpb
    echo -e "PbPb data" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbdata_recoreco_0_20" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbdata_recoreco_20_60" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbdata_recoreco_60_100" >> ${PLOTLIST}_pbpb
    echo -e "hjetshape_final_pbpbdata_recoreco_100_200" >> ${PLOTLIST}_pbpb

    ./plot_results datamc_${1}_${3}_gxi${5}_js_final.root datamc_gxi${5}_${1}_${3}_pbpb ${PLOTLIST}_pbpb 1 $5 $1 $3

    cat ${PLOTLIST}_pp >> $PLOTLIST
    cat ${PLOTLIST}_pbpb >> $PLOTLIST

    ./plot_results datamc_${1}_${3}_gxi${5}_js_final.root datamc_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3

    rm ${PLOTLIST}_pp ${PLOTLIST}_pbpb
elif [ $8 == "rename" ]; then
    hadd -f ${9}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_*_js.root
    ./draw_js ${6} ${9}_${6}_${1}_${3}_gxi${5}_js_merged.root ${9}_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${@:10}

    for i in ${@:10}
    do
      echo -e "$i" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_0_20" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_20_60" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_60_100" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_100_200" >> $PLOTLIST
    done

    ./plot_results ${9}_${6}_${1}_${3}_gxi${5}_js_final.root ${9}_${6}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
else
    hadd -f ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_*_js.root
    ./draw_js ${6} ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${@:8}

    for i in ${@:8}
    do
      echo -e "$i" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_0_20" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_20_60" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_60_100" >> $PLOTLIST
      echo -e "hjetshape_final_${6}_${i}_100_200" >> $PLOTLIST
    done

    ./plot_results ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
fi

rm $PLOTLIST
