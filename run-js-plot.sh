#!/bin/bash

if [ $# -lt 8 ]; then
  echo "Usage: ./run-js-plot.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
  echo "Example: ./run-js-plot.sh 60 1000 30 1 0 pbpbmc closure sgengen sgenreco recogen recoreco"
  exit 1
fi

echo "compiling macros..."
# g++ draw_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_js || exit 1
# g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results || exit 1

PLOTLIST=jsplot_${1}_${3}_${5}_${6}_${7}.list
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
    hadd -f pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pbpb${6}_${1}_${3}_gxi${5}_js_final.root ${7}_pp${6}_${1}_${3}_gxi${5}_js_final.root
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

    ./plot_results pbpbpp${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_pbpb_${9}_pp_${10}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
elif [ $8 = "datamc" ]; then
    hadd -f data_mc_${1}_${3}_gxi${5}_js_merged.root data_ppdata_${1}_${3}_gxi${5}_recoreco_js.root data_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root closure_ppmc_${1}_${3}_gxi${5}_recoreco_js.root closure_pbpbmc_${1}_${3}_gxi${5}_recoreco_js.root
    ./draw_js pbpbdata data_mc_${1}_${3}_gxi${5}_js_merged.root data_mc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppdata data_mc_${1}_${3}_gxi${5}_js_merged.root data_mc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js pbpbmc data_mc_${1}_${3}_gxi${5}_js_merged.root data_mc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppmc data_mc_${1}_${3}_gxi${5}_js_merged.root data_mc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco

    echo -e "pp MC" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_100_200" >> $PLOTLIST
    echo -e "pp data" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb MC" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb data" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results data_mc_${1}_${3}_gxi${5}_js_final.root data_mc_gxi${5}_${1}_${3} $PLOTLIST 0 $5 $1 $3
    rm $PLOTLIST

    echo -e "pp MC" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_ppmc_recoreco_100_200" >> $PLOTLIST
    echo -e "pp data" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_ppdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results data_mc_${1}_${3}_gxi${5}_js_final.root data_mc_gxi${5}_${1}_${3}_pp $PLOTLIST 1 $5 $1 $3
    rm $PLOTLIST

    echo -e "PbPb MC" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbmc_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb data" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hjetshape_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results data_mc_${1}_${3}_gxi${5}_js_final.root data_mc_gxi${5}_${1}_${3}_pbpb $PLOTLIST 1 $5 $1 $3
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

    echo "./plot_results ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3"
    ./plot_results ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${7}_${6}_gxi${5}_${1}_${3} $PLOTLIST 1 $5 $1 $3
fi

rm $PLOTLIST
