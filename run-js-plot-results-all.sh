#!/bin/bash

inputDir="/export/d00/scratch/"$USER"/GJT-out/results/"
inputDirClosure="/export/d00/scratch/"$USER"/GJT-out/results/closure/"

##./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata data $inputDir DUMMY

## sys var plots
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysvar $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysvar $inputDir DUMMY

## sys var plots
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysvar_fitfnc $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysvar_fitfnc $inputDir DUMMY

## sys var plots - All variation together
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysall $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysall $inputDir DUMMY

## sys var plots - All variation together and total systematics
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysalltot $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysalltot $inputDir DUMMY

## sys var plots - All variation together and total systematics in percentage
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysalltotpercnt $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysalltotpercnt $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ratio sysalltotpercnt $inputDir DUMMY

## sys var plots - All photon related variation together and total systematics in percentage
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysallpercnt_phosys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysallpercnt_phosys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ratio sysallpercnt_phosys $inputDir DUMMY

## sys var plots - All jet related variation together and total systematics in percentage
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysallpercnt_jetsys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysallpercnt_jetsys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ratio sysallpercnt_jetsys $inputDir DUMMY

## sys var plots - All remaining variation together and total systematics in percentage
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysallpercnt_othersys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysallpercnt_othersys $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ratio sysallpercnt_othersys $inputDir DUMMY

## closure plots
##./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbmc closure $inputDirClosure DUMMY
##./run-js-plot-results.sh 60 9999 30 1 0 2 ppmc closure $inputDirClosure DUMMY
