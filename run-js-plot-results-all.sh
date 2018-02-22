#!/bin/bash

inputDir="/home/kaya/Documents/EclipseWorkSpace/GJT/gammajettrack/results/"
inputDirClosure=$inputDir"closure/"

##./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata data $inputDir DUMMY

## sys var plots
./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbdata sysvar $inputDir DUMMY
./run-js-plot-results.sh 60 9999 30 1 0 2 ppdata sysvar $inputDir DUMMY

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

## closure plots
##./run-js-plot-results.sh 60 9999 30 1 0 2 pbpbmc closure $inputDirClosure DUMMY
##./run-js-plot-results.sh 60 9999 30 1 0 2 ppmc closure $inputDirClosure DUMMY
