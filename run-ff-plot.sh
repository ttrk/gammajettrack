#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-ff-plot.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [defnFF] [sample] [label] [types...]"
  echo "Example: ./run-ff-plot.sh 60 1000 30 1 0 0 pbpbmc closure sgengen sgenreco recogen recoreco"
  exit 1
fi

echo "phoetmin = $1"
echo "phoetmax = $2"
echo "jetptmin = $3"
echo "trkptmin = $4"
echo "gammaxi  = $5"
echo "defnFF   = $6"
echo "sample   = $7"
echo "label    = $8"
echo "types    = $9"

g++ draw_ff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff.exe || exit 1
echo "g++ draw_ff.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff.exe || exit 1"

g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1
echo "g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1"

PLOTLIST=ffplot_${1}_${3}_${5}_${6}_${7}_${8}.list
if [ -f $PLOTLIST ]; then
    rm $PLOTLIST
fi
touch $PLOTLIST
echo "PLOTLIST : $PLOTLIST"

if [[ $8 == "data" ]]; then
    echo "running $8"
    set -x

    hadd -f data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_ppdata_${1}_${3}_gxi${5}_defnFF${6}_srecoreco_ff.root data_pbpbdata_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root
    ./draw_ff.exe pbpbdata data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root ${1} 0 recoreco
    ./draw_ff.exe ppdata data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root ${1} 0 srecoreco

    echo -e "pp (smeared)" >> $PLOTLIST
    echo -e "hff_final_ppdata_srecoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_ppdata_srecoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_ppdata_srecoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_ppdata_srecoreco_100_200" >> $PLOTLIST
    echo -e "PbPb" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results.exe data_data_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root data_data_gxi${5}_defnFF${6}_${1}_${3} $PLOTLIST 1 $5 $1 $3 2
elif [[ $8 == "datamc" ]]; then
    echo "running $8"
    set -x

    hadd -f data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_ppdata_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root data_pbpbdata_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root closure_ppmc_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root closure_pbpbmc_${1}_${3}_gxi${5}_defnFF${6}_recoreco_ff.root
    ./draw_ff.exe pbpbdata data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root $1 0 recoreco
    ./draw_ff.exe ppdata data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root $1 0 recoreco
    ./draw_ff.exe pbpbmc data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root $1 0 recoreco
    ./draw_ff.exe ppmc data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root $1 0 recoreco

    echo -e "pp MC" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_100_200" >> $PLOTLIST
    echo -e "pp data" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb MC" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb data" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results.exe data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root data_mc_gxi${5}_defnFF${6}_${1}_${3} $PLOTLIST 0 $5 $1 $3 2
    rm $PLOTLIST

    echo -e "pp MC" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_ppmc_recoreco_100_200" >> $PLOTLIST
    echo -e "pp data" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_ppdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results.exe data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root data_mc_gxi${5}_defnFF${6}_${1}_${3}_pp $PLOTLIST 1 $5 $1 $3 2
    rm $PLOTLIST

    echo -e "PbPb MC" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_pbpbmc_recoreco_100_200" >> $PLOTLIST
    echo -e "PbPb data" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_0_20" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_20_60" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_60_100" >> $PLOTLIST
    echo -e "hff_final_pbpbdata_recoreco_100_200" >> $PLOTLIST

    ./plot_results.exe data_mc_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root data_mc_gxi${5}_defnFF${6}_${1}_${3}_pbpb $PLOTLIST 1 $5 $1 $3 2
else
    echo "running $8"
    set -x

    hadd -f ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_*_ff.root
    ./draw_ff.exe ${7} ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_merged.root ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root ${1} 0 ${@:8}

    for i in ${@:8}
    do
      echo -e "$i" >> $PLOTLIST
      echo -e "hff_final_${7}_${i}_0_20" >> $PLOTLIST
      echo -e "hff_final_${7}_${i}_20_60" >> $PLOTLIST
      echo -e "hff_final_${7}_${i}_60_100" >> $PLOTLIST
      echo -e "hff_final_${7}_${i}_100_200" >> $PLOTLIST
    done

    ./plot_results.exe ${8}_${7}_${1}_${3}_gxi${5}_defnFF${6}_ff_final.root ${8}_${7}_gxi${5}_defnFF${6}_${1}_${3} $PLOTLIST 1 $5 $1 $3 2
fi

rm $PLOTLIST
