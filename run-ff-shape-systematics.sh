#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-ff-shape-systematics.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-ff-shape-systematics.sh /export/d00/scratch/$USER/GJT-out/results/sys/ 80 1000 40 1 0 1 pbpbdata jsdata sgengen sgenreco recogen recoreco"
  exit 1
fi

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

if [ $sample = "pbpbdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/PbPb-Data-skim-170911.root"
elif [ $sample = "ppdata" ]; then
    SKIM="/export/d00/scratch/tatar/GJT-out/pp-Data-skim-170911.root"
else
    echo "invalid sample"
    exit 1
fi

echo "compiling macros..."
g++ jetffshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetffshape.exe || exit 1

mkdir -p $outputDir

set -x

if [ $sample = "pbpbmc" ]; then
    SKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_20180115_merged/job0.root"
elif [ $sample = "ppmc" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/pp-MC-skim-180115.root"
elif [ $sample = "pbpbdata" ]; then
    SKIM="/mnt/hadoop/cms/store/user/tatar/GJT-out/skims/PbPb-Data-skim-170911.root"
    MCSKIM="/mnt/hadoop/cms/store/user/katatar/GJT-out/PbPb_MC_skim_20180115_merged/job0.root"
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

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio eta_reflection etagt0p3)

sysIndices="1 2 3 4 6 9 10"
runSysIso=1
if [ $sample = "pbpbdata" ] || [ $sample = "pbpbmc" ]; then
  sysIndices="1 2 3 4 6 9 10 11 12 13 14 15 16"
  runSysIso=1
fi

#set -x

for SYS in ${sysIndices}; do
  echo "running systematics : "${SYSTEMATIC[SYS]}
  for rgLevel in $recogenLevels; do
    #if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
      hiBinMins=(0  20 60  100 0  60)
      hiBinMaxs=(20 60 100 200 60 200)
      if [ $sample = "ppdata" ] && ([ $rgLevel = "recoreco" ] || [ $rgLevel = "corrjsrecoreco" ]); then
        hiBinMins=(100)
        hiBinMaxs=(200)
      fi
      for i1 in ${!hiBinMins[*]}
      do
        hiBinMin=${hiBinMins[i1]}
        hiBinMax=${hiBinMaxs[i1]}
        ./jetffshape.exe $SKIM $sample $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $rgLevel $trkptMin $gammaxi ${outputDir}/${label}_${SYSTEMATIC[SYS]} $SYS $obs &
      done
      wait
    #fi
  done
done
wait

# merge output
for SYS in ${sysIndices}; do
  for rgLevel in $recogenLevels; do
    outputFile=${outputDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root
    inputPrefix=${outputDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}
    hadd -f $outputFile ${inputPrefix}_*_*.root
    rm ${inputPrefix}_*_*.root
  done
  outputFilesMerged=${outputDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
  hadd -f $outputFilesMerged ${outputDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root
done
wait

if [ $runSysIso == 1 ]; then
  echo "running systematics : iso"
  for rgLevel in $recogenLevels; do
    #if [ ! -f ${label}_${sample}_${phoetMin}_${phoetMax}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root ]; then
      hiBinMins=(0  20 60  100 0  60)
      hiBinMaxs=(20 60 100 200 60 200)
      if [ $sample = "ppdata" ] && ([ $rgLevel = "recoreco" ] || [ $rgLevel = "corrjsrecoreco" ]); then
        hiBinMins=(100)
        hiBinMaxs=(200)
      fi
      for i1 in ${!hiBinMins[*]}
      do
        hiBinMin=${hiBinMins[i1]}
        hiBinMax=${hiBinMaxs[i1]}
        ./jetffshape.exe $MCSKIM $MCSAMPLE $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $rgLevel $trkptMin $gammaxi ${outputDir}/${label}_nominal_iso 0 $obs &
        ./jetffshape.exe $MCSKIM $MCSAMPLE $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $rgLevel $trkptMin $gammaxi ${outputDir}/${label}_iso 5 $obs &
      done
      wait
    #fi
  done
  # merge output
  for rgLevel in $recogenLevels; do
    outputFile=${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root
    inputPrefix=${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}
    hadd -f $outputFile ${inputPrefix}_*_*.root
    rm ${inputPrefix}_*_*.root
    outputFile=${outputDir}/${label}_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_${rgLevel}_ffjs.root
    inputPrefix=${outputDir}/${label}_iso_${MCSAMPLE}_${rgLevel}_${phoetMin}_${jetptMin}_${gammaxi}_${obs}
    hadd -f $outputFile ${inputPrefix}_*_*.root
    rm ${inputPrefix}_*_*.root
  done
  outputFilesMerged=${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
  hadd -f $outputFilesMerged ${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root
  outputFilesMerged=${outputDir}/${label}_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
  hadd -f $outputFilesMerged ${outputDir}/${label}_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_*_ffjs.root
fi

