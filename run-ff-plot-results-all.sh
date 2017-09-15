#!/bin/bash

inputDir="/export/d00/scratch/"$USER"/GJT-out/results/"
inputDirClosure="/export/d00/scratch/"$USER"/GJT-out/results/closure/"

./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbdata data $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbdata data $inputDir recoreco
./run-ff-plot-results.sh 80 9999 40 1 0 1 pbpbdata data $inputDir recoreco
./run-ff-plot-results.sh 80 9999 40 1 1 1 pbpbdata data $inputDir recoreco

## sys var plots
./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbdata sysvar $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbdata sysvar $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 0 1 ppdata sysvar $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 ppdata sysvar $inputDir recoreco

#./run-ff-plot-results.sh 80 9999 40 1 0 1 pbpbdata sysvar $inputDir recoreco
#./run-ff-plot-results.sh 80 9999 40 1 1 1 pbpbdata sysvar $inputDir recoreco
#./run-ff-plot-results.sh 80 9999 40 1 0 1 ppdata sysvar $inputDir recoreco
#./run-ff-plot-results.sh 80 9999 40 1 1 1 ppdata sysvar $inputDir recoreco

## sys var plots - All variation together
./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbdata sysall $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbdata sysall $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 0 1 ppdata sysall $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 ppdata sysall $inputDir recoreco

## sys var plots - All variation together and total systematics
./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbdata sysalltot $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbdata sysalltot $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 0 1 ppdata sysalltot $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 ppdata sysalltot $inputDir recoreco

## sys var plots - All variation together and total systematics in percentage
./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbdata sysalltotpercnt $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbdata sysalltotpercnt $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 0 1 ppdata sysalltotpercnt $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 ppdata sysalltotpercnt $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 0 1 ratio sysalltotpercnt $inputDir recoreco
./run-ff-plot-results.sh 60 9999 30 1 1 1 ratio sysalltotpercnt $inputDir recoreco

## closure plots
./run-ff-plot-results.sh 60 9999 30 1 0 1 pbpbmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 60 9999 30 1 1 1 pbpbmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 80 9999 40 1 0 1 pbpbmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 80 9999 40 1 1 1 pbpbmc closure $inputDirClosure DUMMY

./run-ff-plot-results.sh 60 9999 30 1 0 1 ppmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 60 9999 30 1 1 1 ppmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 80 9999 40 1 0 1 ppmc closure $inputDirClosure DUMMY
./run-ff-plot-results.sh 80 9999 40 1 1 1 ppmc closure $inputDirClosure DUMMY
