#!/bin/bash

if [ $# -lt 8 ]; then
    echo "Usage: ./run-js-arithmetic.sh [phoetmin] [phoetmax] [jetptmin] [trkptmin] [gammaxi] [sample] [label] [types...]"
    echo "Example: ./run-js-arithmetic.sh 60 1000 30 1 0 pbpbmc closure sgengen sgenreco recogen recoreco"
    exit 1
fi

echo "compiling macros..."
make draw_js plot_results

set -x

if [ $8 == "data" ]; then
    [ "$9" == "?" ] && echo -e "argument format: ./run-js-arithmetic.sh 60 1000 30 1 0 [pbpb label] [pp label] data [new label]" && exit 1
    [ $# -gt 9 ] && exit 1

    hadd -f ${9}_data_${1}_${3}_gxi${5}_js_merged.root ${7}_ppdata_${1}_${3}_gxi${5}_srecoreco_js.root ${6}_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root
    ./draw_js pbpbdata ${9}_data_${1}_${3}_gxi${5}_js_merged.root ${9}_data_${1}_${3}_gxi${5}_js_final.root ${1} 0 recoreco
    ./draw_js ppdata ${9}_data_${1}_${3}_gxi${5}_js_merged.root ${9}_data_${1}_${3}_gxi${5}_js_final.root ${1} 0 srecoreco

    ./quick-js-plot.sh ${@:1:5} dummy $9 final data
elif [ $8 == "datamc" ]; then
    [ "$9" == "?" ] && echo -e "argument format: ./run-js-arithmetic.sh 60 1000 30 1 0 [data label] [mc label] datamc [new label]" && exit 1
    [ $# -gt 9 ] && exit 1

    hadd -f ${9}_datamc_${1}_${3}_gxi${5}_js_merged.root ${6}_ppdata_${1}_${3}_gxi${5}_recoreco_js.root ${6}_pbpbdata_${1}_${3}_gxi${5}_recoreco_js.root ${7}_ppmc_${1}_${3}_gxi${5}_recoreco_js.root ${7}_pbpbmc_${1}_${3}_gxi${5}_recoreco_js.root
    ./draw_js pbpbdata ${9}_datamc_${1}_${3}_gxi${5}_js_merged.root ${9}_datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppdata ${9}_datamc_${1}_${3}_gxi${5}_js_merged.root ${9}_datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js pbpbmc ${9}_datamc_${1}_${3}_gxi${5}_js_merged.root ${9}_datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco
    ./draw_js ppmc ${9}_datamc_${1}_${3}_gxi${5}_js_merged.root ${9}_datamc_${1}_${3}_gxi${5}_js_final.root $1 0 recoreco

    ./quick-js-plot.sh ${@:1:5} dummy $9 final datamc
elif [ $8 == "rename" ]; then
    [ "$9" == "?" ] && echo -e "argument format: ./run-js-arithmetic.sh 60 1000 30 1 0 [sample] [old label] rename [new label]" && exit 1

    [ -f ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ] &&
        cp ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ${9}_${6}_${1}_${3}_gxi${5}_js_merged.root ||
        hadd -f ${9}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_*_js.root
    ./draw_js ${6} ${9}_${6}_${1}_${3}_gxi${5}_js_merged.root ${9}_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${@:10}

    ./quick-js-plot.sh ${@:1:6} $9 final ${@:10}
else
    hadd -f ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_*_js.root
    ./draw_js ${6} ${7}_${6}_${1}_${3}_gxi${5}_js_merged.root ${7}_${6}_${1}_${3}_gxi${5}_js_final.root ${1} 0 ${@:8}

    ./quick-js-plot.sh ${@:1:7} final ${@:8}
fi
