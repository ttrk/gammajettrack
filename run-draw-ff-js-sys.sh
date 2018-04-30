#!/bin/bash

if [ $# -lt 9 ]; then
  echo "Usage: ./run-draw-ff-js-sys.sh [outputDir] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [obs] [sample] [label] [types...]"
  echo "Example: ./run-draw-ff-js-sys.sh /export/d00/scratch/"$USER"/GJT-out/results/sys/ 60 1000 30 1 0 2 pbpbdata jssys recoreco"
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

g++ draw_ff_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff_js.exe || exit 1
echo "g++ draw_ff_js.C $(root-config --cflags --libs) -Werror -Wall -O2 -o draw_ff_js.exe || exit 1"

if [ $sample = "pbpbdata" ]; then
    MCSAMPLE="pbpbmc"
elif [ $sample = "ppdata" ]; then
    MCSAMPLE="ppmc"
else
    echo "invalid sample"
    exit 1
fi

SYSTEMATIC=(placeholder jes_up jes_down jer pes iso ele_rej purity_up purity_down tracking_up tracking_down jes_qg_up jes_qg_down longrange tracking_ratio eta_reflection etagt0p3 phoeffcorr)

sysIndices="1 2 3 4 6 9 10 13 17"
runSysIso=1
if [ $sample = "pbpbdata" ] || [ $sample = "pbpbmc" ]; then
  sysIndices="1 2 3 4 6 9 10 12 13 14 15 16 17"
  runSysIso=1
fi

purityGroup=0

set -x
for SYS in ${sysIndices}; do
  echo "running systematics : " ${SYSTEMATIC[SYS]}
  inputFile=${outputDir}/${label}_${SYSTEMATIC[SYS]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
  if [ -f $inputFile ]; then
    outputFile="${inputFile/_merged.root/_final.root}"
    ./draw_ff_js.exe $sample $inputFile $outputFile $phoetMin $purityGroup $recogenLevels
  else
    echo "inputFile : $inputFile does not exist. Skipping this systematics."
  fi
done
wait

## special systematics
# purity
indexPurityUp=7
indexPurityDown=8
# construct the input directory for nominal results
inputDirNominal="${outputDir/\/sys/}"	# remove "/sys" from dirname
inputFile=${inputDirNominal}/jsdata_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
if [ -f $inputFile ]; then
  outputFilePurityUp=${outputDir}/${label}_${SYSTEMATIC[indexPurityUp]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  outputFilePurityDown=${outputDir}/${label}_${SYSTEMATIC[indexPurityDown]}_${sample}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  ./draw_ff_js.exe $sample $inputFile $outputFilePurityUp $phoetMin 2 $recogenLevels
  ./draw_ff_js.exe $sample $inputFile $outputFilePurityDown $phoetMin -2 $recogenLevels
else
  echo "inputFile : $inputFile does not exist. Skipping this systematics."
fi

# isolation

# construct the input directory for nominal results
indexIso=5
inputFileNominalIso=${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
inputFileIso=${outputDir}/${label}_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_merged.root
if [ -f $inputFileNominalIso ] && [ -f $inputFileIso ]; then
  outputFileNominalIso=${outputDir}/${label}_nominal_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  outputFileIso=${outputDir}/${label}_iso_${MCSAMPLE}_${phoetMin}_${jetptMin}_gxi${gammaxi}_obs${obs}_ffjs_final.root
  ./draw_ff_js.exe $MCSAMPLE $inputFileNominalIso $outputFileNominalIso $phoetMin $purityGroup $recogenLevels
  ./draw_ff_js.exe $MCSAMPLE $inputFileIso $outputFileIso $phoetMin $purityGroup $recogenLevels
else
  echo "inputFile : $inputFileNominalIso and/or $inputFileIso does not exist. Skipping this systematics."
fi

