#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: ./run-calc-spectra-weights-ffjs.sh [input] [output]"
  echo "Example: ./run-calc-spectra-weights-ffjs.sh /export/d00/scratch/tatar/GJT-out/results//jsdata_data_60_30_gxi0_obs2_ffjs_final.root /export/d00/scratch/tatar/GJT-out/results//jsdata_data_60_30_gxi0_obs2_ffjs_spectra_weights.root"
  exit 1
fi

echo "input = $1"
echo "output = $2"

inputpp=$1
inputpbpb=$1
output=$2

g++ calc_spectra_weights.C $(root-config --cflags --libs) -Werror -Wall -O2 -o calc_spectra_weights.exe || exit 1

set -x

./calc_spectra_weights.exe $inputpp $inputpbpb $output

