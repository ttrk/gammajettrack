#!/bin/bash

if [ $# -lt 10 ]; then
  echo "Usage: ./run-js-closure-condor.sh [input] [output] [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
  echo "Example: ./run-js-closure-condor.sh [input] [output] 80 1000 40 1 0 pbpbmc jsclosure sgengen sgenreco recogen recoreco"
  exit 1
fi

echo "compiling macros..."
g++ jetshape.C $(root-config --cflags --libs) -Werror -Wall -O2 -o jetshape || exit 1

NJOBS=$(ls -l ${1}/*.root | wc -l)

CONDOR_INPUT=inputs_for_condor_${3}_${5}_${7}_${8}_${9}
mkdir $CONDOR_INPUT
cp jetshape $CONDOR_INPUT
cp PbPb-weights.root $CONDOR_INPUT
cp pp-weights.root $CONDOR_INPUT

set -x

./run-script-on-condor.sh $NJOBS js-closure-condor.sh $CONDOR_INPUT $2 $1 ${@:3}
