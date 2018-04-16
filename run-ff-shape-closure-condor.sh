#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-ff-shape-closure-condor.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-ff-shape-closure-condor.sh /mnt/hadoop/cms/store/user/${USER}/GJT-out/results/closure/ 80 1000 40 1 0 1 pbpbmc ffclosure sgengen sgenreco recogen recoreco &> logFile &"
  exit 1
fi

# transfer the standard out to a log file and then grep the condor submission commands via
# grep "condor_submit" <logfile>

outputDir=$1
phoetMin=$2
phoetMax=$3
jetptMin=$4
trkptMin=$5
gammaxi=$6
obs=$7
sample=$8
label=$9
recogenLevels=${@:10}

echo "outputDir = $outputDir"
echo "phoetMin  = $phoetMin"
echo "phoetMax  = $phoetMax"
echo "jetptMin  = $jetptMin"
echo "trkptMin  = $trkptMin"
echo "gammaxi   = $gammaxi"
echo "obs       = $obs"
echo "sample    = $sample"
echo "label     = $label"
echo "recogenLevels = $recogenLevels"

if [[ $outputDir != /mnt/hadoop/* ]]; then
    echo "output directory must be under /mnt/hadoop/"
    exit 1
fi

if [ $sample = "pbpbmc" ]; then
    SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_20180115_merged/job0.root"
    if [[ $outputDir = *flt50ext* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_Flt50_ext_skim_20180413_merged/job0.root"
    elif [[ $outputDir = *flt50* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_Flt50_skim_20180413_merged/job0.root"
    elif [[ $outputDir = *0302* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_20180302_merged/job0.root"
    elif [[ $outputDir = *0309* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_20180309_merged/job0.root"
    elif [[ $outputDir = *photonOnly* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_photonOnly_merged/job0.root"
    fi
elif [ $sample = "ppmc" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/pp-MC-skim-180115.root"
    if [[ $outputDir = *stdPho* ]]; then
      SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/pp_MC_stdPho_skim_20180309_merged/job0.root"
    fi
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1

#set -x

echo prepare condor jobs for closure histograms
for rgLevel in $recogenLevels; do
  if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
    hiBinMins=(0  20 60  100)
    hiBinMaxs=(20 60 100 200)
    if [ $sample = "ppmc" ]; then
      hiBinMins=(100)
      hiBinMaxs=(200)
    fi
    for i1 in ${!hiBinMins[*]}
    do
      hiBinMin=${hiBinMins[i1]}
      hiBinMax=${hiBinMaxs[i1]}
     ./condor/condorSubmit_jetffshape.sh $SKIM $sample $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $trkptMin $gammaxi $label 0 $obs $rgLevel $outputDir
    done
  fi
done
wait

#for rgLevel in $recogenLevels; do
#  if [ ! -f ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
#    hadd -f ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ${label}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
#    rm ${label}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}_*_*.root
#  fi
#done

#./run-ff-shape-plot.sh $@

#outDir="/export/d00/scratch/"$USER"/GJT-out/results/closure/"
#mkdir -p $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root $outDir
#mv ${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root $outDir
