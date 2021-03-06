#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-js-plot-results.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [inputDir] [types...]"
  echo "Example: ./run-js-plot-results.sh 60 1000 30 1 0 0 pbpbmc closure results/closure/ sgengen sgenreco recogen recoreco"
  echo "Example: ./run-js-plot-results.sh 60 9999 30 1 0 1 pbpbdata data results/ recoreco"
  exit 1
fi

phoetmin=$1
phoetmax=$2
jetptmin=$3
trkptmin=$4
gammaxi=$5
obs=$6
sample=$7
label=$8
inputDir=$9
types=${10}

echo "phoetmin = $phoetmin"
echo "phoetmax = $phoetmax"
echo "jetptmin = $phoetmin"
echo "trkptmin = $trkptmin"
echo "gammaxi  = $gammaxi"
echo "obs      = $obs"
echo "sample   = $sample"
echo "label    = $label"
echo "inputDir = $inputDir"
echo "types    = ${types}"

g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1
echo "g++ plot_results.C $(root-config --cflags --libs) -Werror -Wall -O2 -o plot_results.exe || exit 1"

inDir=$inputDir
outDir=$inDir
mkdir -p $outDir
echo "inDir : $inDir"
echo "outDir : $outDir"
plotOption=1

PLOTLIST=jsplot_${phoetmin}_${jetptmin}_${gammaxi}_${obs}_${sample}_${label}.list
if [ -f $PLOTLIST ]; then
    rm $PLOTLIST
fi
touch $PLOTLIST
echo "PLOTLIST : $PLOTLIST"

if [[ $label == "data" ]]; then
    echo "running $label"
    set -x
    inPrefix=$inDir"jsdata"
    outPrefix=$outDir"jsdata"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    echo -e "pp" >> $PLOTLIST
    echo -e "hjs_final_ppdata_corrjsrecoreco_100_200" >> $PLOTLIST
    echo -e "hjs_final_ppdata_corrjsrecoreco_100_200" >> $PLOTLIST
    echo -e "hjs_final_ppdata_corrjsrecoreco_100_200" >> $PLOTLIST
    echo -e "hjs_final_ppdata_corrjsrecoreco_100_200" >> $PLOTLIST
    echo -e "PbPb" >> $PLOTLIST
    echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200" >> $PLOTLIST
    echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100" >> $PLOTLIST
    echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60" >> $PLOTLIST
    echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20" >> $PLOTLIST

    nCols=4
    ./plot_results.exe ${inPrefix}_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_data_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption $sysFile $nCols
elif [[ $label == "closure" ]]; then
    echo "running $label"
    set -x
    inPrefix=$inDir"jsclosure"
    outPrefix=$outDir"jsclosure"

    if [[ $sample == "pbpbmc" ]]; then
      echo -e "smeared gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_0_20" >> $PLOTLIST
      echo -e "reco jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_recoreco_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_recoreco_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_recoreco_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_recoreco_0_20" >> $PLOTLIST

      nCols=4
      ./plot_results.exe ${inPrefix}_pbpbmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_pbpbmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST

      echo -e "smeared signal gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgen0gen0_0_20" >> $PLOTLIST
      echo -e "smeared gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgengen_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgengen_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgengen_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_sgengen_0_20" >> $PLOTLIST

      nCols=4
      ./plot_results.exe ${inPrefix}_pbpbmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_pbpbmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin}_sgengen $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST

      echo -e "signal gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gen0gen0_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gen0gen0_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gen0gen0_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gen0gen0_0_20" >> $PLOTLIST
      echo -e "gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gengen_100_200" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gengen_60_100" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gengen_20_60" >> $PLOTLIST
      echo -e "hjs_final_pbpbmc_gengen_0_20" >> $PLOTLIST

      nCols=4
      ./plot_results.exe ${inPrefix}_pbpbmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_pbpbmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin}_gengen $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ppmc" ]]; then
      echo -e "gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_100_200" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_60_100" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_20_60" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_0_20" >> $PLOTLIST
      echo -e "reco jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_100_200" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_60_100" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_20_60" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_0_20" >> $PLOTLIST

      nCols=4
      ./plot_results.exe ${inPrefix}_ppmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_ppmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "closure4" ]]; then
    echo "running $label"
    set -x
    inPrefix=$inDir"jsclosure"
    outPrefix=$outDir"jsclosure"

    if [[ $sample == "pbpbmc" ]]; then
      echo -e "sgengen" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgengen_100_200" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgengen_60_100" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgengen_20_60" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgengen_0_20" >> $PLOTLIST 
      echo -e "sgenreco" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgenreco_100_200" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgenreco_60_100" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgenreco_20_60" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_sgenreco_0_20" >> $PLOTLIST 
      echo -e "recogen" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recogen_100_200" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recogen_60_100" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recogen_20_60" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recogen_0_20" >> $PLOTLIST 
      echo -e "recoreco" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recoreco_100_200" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recoreco_60_100" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recoreco_20_60" >> $PLOTLIST 
      echo -e "hjs_final_pbpbmc_recoreco_0_20" >> $PLOTLIST 

      nCols=4
      ./plot_results.exe ${inPrefix}_pbpbmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_pbpbmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin}_4 $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ppmc" ]]; then
      echo -e "Pythia gen jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_100_200" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_60_100" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_20_60" >> $PLOTLIST
      echo -e "hjs_final_ppmc_sgengen_0_20" >> $PLOTLIST
      echo -e "Pythia reco jets, tracks" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_100_200" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_60_100" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_20_60" >> $PLOTLIST
      echo -e "hjs_final_ppmc_recoreco_0_20" >> $PLOTLIST

      nCols=4
      ./plot_results.exe ${inPrefix}_ppmc_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${outPrefix}_ppmc_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "sysvar" ]]; then
    echo "running $label"
    set -x
    inPrefix=$inDir"data_sysvar"
    outPrefix=$outDir"data_sysvar"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    sysVarIndices=(0 2 3 4 5 6 8 10 11 12)
    sysVarLabels=(jes_up      jes_down jer pes iso ele_rej purity_up        purity_down  tracking_up      tracking_down jes_qg_down jes_ue noUEscale)
    sysVarNames=( jes_up_down NULL     jer pes iso ele_rej purity_up_down   NULL         tracking_up_down NULL          jes_qg_down jes_ue noUEscale)
    sysVarUpDowns=(1 1 0 0 0 0 1 1 1 1 0 0 0)
    sysVarTitles=(
    "JES UP"
    "JES DOWN"
    "JER"
    "PES"
    "isolation"
    "electron rejection"
    "purity UP"
    "purity DOWN"
    "tracking UP"
    "tracking DOWN"
    "JES quark-jet"
    "JES data-MC UE diff"
    "bkg sub, UE scaling"
     )

    for i1 in ${!sysVarIndices[*]}
    do
      iSys=${sysVarIndices[i1]}
      sysVarLabel=${sysVarLabels[iSys]}
      sysVarName=${sysVarNames[iSys]}
      sysVarUpDown=${sysVarUpDowns[iSys]}
      sysVarTitle=${sysVarTitles[iSys]}
      if [[ $sample == "pbpbdata" ]]; then
        echo -e "PbPb - nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_variation" >> $PLOTLIST
        if [[ $sysVarUpDown > 0 ]]; then
          echo -e ${sysVarTitles[$((iSys+1))]} >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
        fi

        nCols=4
        ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_${sysVarName}_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      elif [[ $sample == "ppdata" ]]; then
        echo -e "pp - nominal" >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        if [[ $sysVarUpDown > 0 ]]; then
          echo -e ${sysVarTitles[$((iSys+1))]} >> $PLOTLIST
          echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
        fi

        nCols=1
        ./plot_results.exe $sysFile ${outPrefix}_ppdata_${sysVarName}_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      fi
      rm $PLOTLIST
    done
elif [[ $label == "sysvar_fitfnc" ]]; then
    echo "running $label"
    set -x
    inPrefix=$inDir"data_sysvar_fitfnc"
    outPrefix=$outDir"data_sysvar_fitfnc"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    sysVarIndices=(0 2 3 4 5 6 8 10 11 12 13)
    sysVarLabels=(jes_up      jes_down jer pes iso ele_rej purity_up        purity_down  tracking_up      tracking_down jes_qg_down phoeffcorr jes_ue noUEscale)
    sysVarNames=( jes_up_down NULL     jer pes iso ele_rej purity_up_down   NULL         tracking_up_down NULL          jes_qg_down phoeffcorr jes_ue noUEscale)
    sysVarUpDowns=(1 1 0 0 0 0 1 1 1 1 0 0 0 0)
    sysVarTitles=(
    "JES UP"
    "JES DOWN"
    "JER"
    "PES"
    "isolation"
    "electron rejection"
    "purity UP"
    "purity DOWN"
    "tracking UP"
    "tracking DOWN"
    "JES quark-jet"
    "photon efficiency"
    "JES data-MC UE diff"
    "bkg sub, UE scaling"
     )

    for i1 in ${!sysVarIndices[*]}
    do
      iSys=${sysVarIndices[i1]}
      sysVarLabel=${sysVarLabels[iSys]}
      sysVarName=${sysVarNames[iSys]}
      sysVarUpDown=${sysVarUpDowns[iSys]}
      sysVarTitle=${sysVarTitles[iSys]}
      if [[ $sample == "pbpbdata" ]]; then
        echo -e "PbPb - nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_variation" >> $PLOTLIST
        if [[ $sysVarUpDown > 0 ]]; then
          echo -e ${sysVarTitles[$((iSys+1))]} >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
          echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
        fi

        nCols=4
        ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_${sysVarName}_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      elif [[ $sample == "ppdata" ]]; then
        echo -e "pp - nominal" >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_nominal" >> $PLOTLIST
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        if [[ $sysVarUpDown > 0 ]]; then
          echo -e ${sysVarTitles[$((iSys+1))]} >> $PLOTLIST
          echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"${sysVarLabels[$((iSys+1))]}"_variation" >> $PLOTLIST
        fi

        nCols=1
        ./plot_results.exe $sysFile ${outPrefix}_ppdata_${sysVarName}_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      fi
      rm $PLOTLIST
    done
elif [[ $label == "sysall" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      #sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10)
      #sysVarLabels=(jes_qg_down jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange js_nonclosure)
      sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10 11)
      sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down js_nonclosure jes_qg_down jes_ue noUEscale)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "PES"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "nonclosure"
      "JES quark-jet"
      "JES data-MC UE diff"
      "bkg sub, UE scaling"
       )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_variation" >> $PLOTLIST
      done

    nCols=4
    ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysall_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
    rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      #sysVarIndices=(0 1 2 3 4 5 6 7 8)
      #sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange)
      sysVarIndices=(0 1 2 3 4 5 6 7 8)
      sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down js_nonclosure)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "PES"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "js_nonclosure"
       )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
      done

    nCols=1
    ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysall_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
    rm $PLOTLIST
    fi
elif [[ $label == "sysalltot" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      #sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10)
      #sysVarLabels=(jes_qg_down jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange js_nonclosure)
      sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10 11)
      sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down js_nonclosure jes_qg_down jes_ue noUEscale)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "PES"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "nonclosure"
      "JES quark-jet"
      "JES data-MC UE diff"
      "bkg sub, UE scaling"
       )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_variation" >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_variation" >> $PLOTLIST
      done
      echo -e "total sys" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_systematics" >> $PLOTLIST

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysalltot_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      #sysVarIndices=(0 1 2 3 4 5 6 7 8)
      #sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange)
      sysVarIndices=(0 1 2 3 4 5 6 7 8)
      sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down js_nonclosure)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "PES"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "nonclosure"
       )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_variation" >> $PLOTLIST
      done
      echo -e "total sys" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_systematics" >> $PLOTLIST

      nCols=1
      ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysalltot_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "sysalltotpercnt" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      #sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10)
      #sysVarLabels=(jes_qg_down jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange js_nonclosure)
      sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
      sysVarLabels=(jes_down jes_up jer iso ele_rej purity_up purity_down tracking_up tracking_down js_nonclosure pes phoeffcorr jes_qg_down bkgsub jes_ue noUEscale)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "tracking UP"
      "tracking DOWN"
      "nonclosure"
      "PES"
      "photon efficiency"
      "JES quark-jet"
      "bkg sub, #eta-reflection"
      "JES data-MC UE diff"
      "bkg sub, UE scaling"
       )
      sysMethodIndices=(1 1 1 1 1 1 1 1 1 0 1 1 1 0 0 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done
      echo -e "total sys" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_systematics" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_systematics" >> $PLOTLIST

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysalltotpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10)
      #sysVarLabels=(jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange)
      sysVarLabels=(jes_down jes_up jer iso ele_rej purity_up purity_down tracking_up tracking_down js_nonclosure phoeffcorr)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "tracking UP"
      "tracking DOWN"
      "nonclosure"
      "photon efficiency"
       )
      sysMethodIndices=(1 1 1 1 1 1 1 1 1 0 1)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done
      echo -e "total sys" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_systematics" >> $PLOTLIST

      nCols=1
      ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysalltotpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ratio" ]]; then

      #sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10)
      #sysVarLabels=(jes_qg_down jes_down jes_up jer pes iso ele_rej purity_up purity_down longrange js_nonclosure)
      sysVarIndices=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14)
      sysVarLabels=(jes_down jes_up jer iso ele_rej purity_up purity_down tracking_ratio js_nonclosure pes phoeffcorr jes_qg_down bkgsub jes_ue noUEscale)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      #"Long range"
      "tracking eff."
      "nonclosure"
      "PES"
      "photon efficiency"
      "JES quark-jet"
      "bkg sub, #eta-reflection"
      "JES data-MC UE diff"
      "bkg sub, UE scaling"
       )
      sysMethodIndices=(1 1 1 1 1 1 1 0 0 1 1 1 0 0 0)
      sysMethodSuffices=(
      "ratio"
      "ratio_fit"
      )

      echo -e "PbPb/pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ratio_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done
      echo -e "total sys" >> $PLOTLIST
      echo -e "hjs_final_ratio_100_200_systematics" >> $PLOTLIST
      echo -e "hjs_final_ratio_60_100_systematics" >> $PLOTLIST
      echo -e "hjs_final_ratio_20_60_systematics" >> $PLOTLIST
      echo -e "hjs_final_ratio_0_20_systematics" >> $PLOTLIST

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_ratio_sysalltotpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "sysallpercnt_phosys" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar_phosys"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      sysVarIndices=(0 1 2 3 4 5)
      sysVarLabels=(iso ele_rej purity_up purity_down pes phoeffcorr)
      sysVarTitles=(
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      "PES"
      "photon efficiency"
       )
      sysMethodIndices=(1 1 1 1 1 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      sysVarIndices=(0 1 2 3)
      sysVarLabels=(iso ele_rej purity_up purity_down)
      sysVarTitles=(
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
       )
      sysMethodIndices=(1 1 1 1)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=1
      ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ratio" ]]; then

      sysVarIndices=(0 1 2 3 4 5)
      sysVarLabels=(iso ele_rej purity_up purity_down pes phoeffcorr)
      sysVarTitles=(
      "isolation"
      "electron rejection"
      "purity UP"
      "purity DOWN"
      "PES"
      "photon efficiency"
       )
      sysMethodIndices=(1 1 1 1 1 0)
      sysMethodSuffices=(
      "ratio"
      "ratio_fit"
      )

      echo -e "PbPb/pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ratio_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_ratio_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "sysallpercnt_jetsys" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar_jetsys"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      sysVarIndices=(0 1 2 3 4)
      sysVarLabels=(jes_down jes_up jer jes_qg_down jes_ue)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "JES quark-jet"
      "JES data-MC UE diff"
       )
      sysMethodIndices=(1 1 1 1 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      sysVarIndices=(0 1 2)
      sysVarLabels=(jes_down jes_up jer)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
       )
      sysMethodIndices=(1 1 1)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=1
      ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ratio" ]]; then
      sysVarIndices=(0 1 2 3 4)
      sysVarLabels=(jes_down jes_up jer jes_qg_down jes_ue)
      sysVarTitles=(
      "JES DOWN"
      "JES UP"
      "JER"
      "JES quark-jet"
      "JES data-MC UE diff"
       )
      sysMethodIndices=(1 1 1 1 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb/pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ratio_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_ratio_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
elif [[ $label == "sysallpercnt_othersys" ]]; then
    echo "running $label"
    set -x
    outPrefix=$outDir"data_sysvar_othersys"

    sysFile=${inDir}sys/jssys_data_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs-systematics.root
    echo "sysFile : $sysFile"

    if [[ $sample == "pbpbdata" ]]; then
      sysVarIndices=(0 1 2 3 4)
      sysVarLabels=(js_nonclosure tracking_up tracking_down bkgsub noUEscale)
      sysVarTitles=(
      "nonclosure"
      "tracking UP"
      "tracking DOWN"
      "bkg sub, #eta-reflection"
      "bkg sub, UE scaling"
       )
      sysMethodIndices=(0 1 1 0 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb - nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_pbpbdata_corrjsrecoreco_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_pbpbdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    
    elif [[ $sample == "ppdata" ]]; then
      sysVarIndices=(0 1 2)
      sysVarLabels=(js_nonclosure tracking_up tracking_down)
      sysVarTitles=(
      "nonclosure"
      "tracking UP"
      "tracking DOWN"
       )
      sysMethodIndices=(0 1 1)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ppdata_corrjsrecoreco_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=1
      ./plot_results.exe $sysFile ${outPrefix}_ppdata_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    elif [[ $sample == "ratio" ]]; then
      sysVarIndices=(0 1 2 3)
      sysVarLabels=(js_nonclosure tracking_ratio bkgsub noUEscale)
      sysVarTitles=(
      "nonclosure"
      "tracking eff."
      "bkg sub, #eta-reflection"
      "bkg sub, UE scaling"
       )
      sysMethodIndices=(0 1 0 0)
      sysMethodSuffices=(
      "ratio" 
      "ratio_fit"
      )

      echo -e "PbPb/pp - nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_100_200_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_60_100_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_20_60_nominal" >> $PLOTLIST
      echo -e "hjs_final_ratio_0_20_nominal" >> $PLOTLIST
      for i1 in ${!sysVarIndices[*]}
      do
        iSys=${sysVarIndices[i1]}
        sysVarLabel=${sysVarLabels[iSys]}
        sysVarTitle=${sysVarTitles[iSys]}
        iSysMethod=${sysMethodIndices[i1]}
        sysMethodSuffix=${sysMethodSuffices[iSysMethod]}
        echo -e $sysVarTitle >> $PLOTLIST
        echo -e "hjs_final_ratio_100_200_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_60_100_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_20_60_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
        echo -e "hjs_final_ratio_0_20_"$sysVarLabel"_"$sysMethodSuffix >> $PLOTLIST
      done

      nCols=4
      ./plot_results.exe $sysFile ${outPrefix}_ratio_sysallpercnt_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
      rm $PLOTLIST
    fi
else
    echo "running $label"
    set -x

    for i in ${@:10}
    do
      echo -e "$i" >> $PLOTLIST
      echo -e "hjs_final_${sample}_${i}_0_20" >> $PLOTLIST
      echo -e "hjs_final_${sample}_${i}_20_60" >> $PLOTLIST
      echo -e "hjs_final_${sample}_${i}_60_100" >> $PLOTLIST
      echo -e "hjs_final_${sample}_${i}_100_200" >> $PLOTLIST
    done

    nCols=4
    ./plot_results.exe ${label}_${sample}_${phoetmin}_${jetptmin}_gxi${gammaxi}_obs${obs}_ffjs_final.root ${label}_${sample}_gxi${gammaxi}_obs${obs}_${phoetmin}_${jetptmin} $PLOTLIST 1 $gammaxi $phoetmin $phoetmin $plotOption DUMMYSYS $nCols
fi

rm $PLOTLIST
