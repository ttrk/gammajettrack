#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-ff-shape-systematics-condor.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-ff-shape-systematics-condor.sh /mnt/hadoop/cms/store/user/${USER}/GJT-out/results/closure/ 80 1000 40 1 0 1 pbpbmc ffclosure sgengen sgenreco recogen recoreco &> logFile &"
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
    SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_Flt50_ext_skim_20180413_merged/job0.root"
elif [ $sample = "ppmc" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/pp-MC-skim-180115.root"
elif [ $sample = "pbpbdata" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/skims/PbPb-Data-skim-170911.root"
    MCSKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_Flt50_ext_skim_20180413_merged/job0.root"
    MCSAMPLE="pbpbmc"
elif [ $sample = "ppdata" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/skims/pp-Data-skim-170911.root"
    MCSKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/pp-MC-skim-180115.root"
    MCSAMPLE="ppmc"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio eta_reflection etagt0p3 phoeffcorr jes_ue_up jes_ue_down noUEscale)

sysIndices="1 2 3 4 6 9 10 17"
runSysIso=1
if [ $sample = "pbpbdata" ] || [ $sample = "pbpbmc" ]; then
  sysIndices="1 2 3 4 6 9 10 12 13 14 15 16 17 18 19 20"
  runSysIso=1
fi

#set -x

echo "prepare condor jobs"
for SYS in ${sysIndices}; do
  tmpSKIM=$SKIM
  tmpSample=$sample
  if [ $SYS == 18 ] || [ $SYS == 19 ]; then
    tmpSKIM=$MCSKIM
    tmpSample=$MCSAMPLE
  fi
  for rgLevel in $recogenLevels; do
    #if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
      hiBinMins=(0  20 60  100 0  60)
      hiBinMaxs=(20 60 100 200 60 200)
      for i1 in ${!hiBinMins[*]}
      do
        hiBinMin=${hiBinMins[i1]}
        hiBinMax=${hiBinMaxs[i1]}
       ./condor/condorSubmit_jetffshape.sh $tmpSKIM $tmpSample $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $trkptMin $gammaxi ${label}_${SYSTEMATIC[SYS]} $SYS $obs $rgLevel $outputDir
      done
    #fi
  done
done
wait

if [ $runSysIso == 1 ]; then
  for rgLevel in $recogenLevels; do
    #if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
      hiBinMins=(0  20 60  100 0  60)
      hiBinMaxs=(20 60 100 200 60 200)
      for i1 in ${!hiBinMins[*]}
      do
        hiBinMin=${hiBinMins[i1]}
        hiBinMax=${hiBinMaxs[i1]}
       ./condor/condorSubmit_jetffshape.sh $MCSKIM $MCSAMPLE $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $trkptMin $gammaxi ${label}_nominal_iso 0 $obs $rgLevel $outputDir
       ./condor/condorSubmit_jetffshape.sh $MCSKIM $MCSAMPLE $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $trkptMin $gammaxi ${label}_iso 5 $obs $rgLevel $outputDir
      done
    #fi
  done
fi


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
