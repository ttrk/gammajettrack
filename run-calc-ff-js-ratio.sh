#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-calc-ff-js-ratio.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-calc-ff-js-ratio.sh /export/d00/scratch/"$USER"/GJT-out/results/ 60 1000 30 1 0 2 data jsdata sgengen sgenreco recogen recoreco"
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

mkdir -p $outputDir

g++ calc_ff_js_ratio.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ff_js_ratio.exe || exit 1
echo "g++ calc_ff_js_ratio.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_ff_js_ratio.exe || exit 1"

set -x

finalAll=${outputDir}/${label}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
finalPP=${outputDir}/${label}_pp${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
finalPBPB=${outputDir}/${label}_pbpb${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root

hadd -f $finalAll $finalPP $finalPBPB
for rgLevel in $recogenLevels; do
  rgLevelPP=$rgLevel
  if [ $rgLevelPP = "recoreco" ]; then
    rgLevelPP="srecoreco"
  fi
  ./calc_ff_js_ratio.exe $finalAll pp${sample}_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll
  ./calc_ff_js_ratio.exe $finalAll pp${sample}_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll hjs_normJet
done
wait

finalAllreweight=${outputDir}/${label}_${sample}reweight_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
finalPPreweight=${outputDir}/${label}_pp${sample}reweight_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root

hadd -f $finalAllreweight $finalPPreweight $finalPBPB
for rgLevel in $recogenLevels; do
  rgLevelPP=$rgLevel
  ./calc_ff_js_ratio.exe $finalAll pp${sample}reweight_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll
  ./calc_ff_js_ratio.exe $finalAll pp${sample}reweight_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll hjs_normJet
done
wait

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrangecalc tracking_ratio eta_reflection etagt0p3)
SYSTEMATICPP=(NOM       jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down NOM       NOM         longrangecalc NOM            NOM            NOM)

sysIndices="1 2 3 4 6 9 10 11 12 13 14 15 16"
runSysIso=1

for SYS in ${sysIndices}; do
  syslabel=${SYSTEMATIC[SYS]}
  syslabelPP=${SYSTEMATICPP[SYS]}
  echo "running ratio for systematics : "${syslabel}
  label="jssys"
  finalAll=${outputDir}/sys/${label}_${syslabel}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  finalPP=${outputDir}/sys/${label}_${syslabelPP}_pp${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  if [ $syslabelPP = NOM ]; then
    finalPP=${outputDir}/jsdata_pp${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  fi 
  finalPBPB=${outputDir}/sys/${label}_${syslabel}_pbpb${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  hadd -f $finalAll $finalPP $finalPBPB
  for rgLevel in $recogenLevels; do
    rgLevelPP=$rgLevel
    if [ $rgLevelPP = "recoreco" ]; then
      rgLevelPP="srecoreco"
    fi
    ./calc_ff_js_ratio.exe $finalAll pp${sample}_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll
    ./calc_ff_js_ratio.exe $finalAll pp${sample}_${rgLevelPP} pbpb${sample}_${rgLevel} $finalAll hjs_normJet
  done
done
wait

