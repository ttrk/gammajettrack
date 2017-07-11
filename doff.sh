
if [ $# -lt 6 ]
then
  echo "Usage: ./doff.sh <phoetmin> <phoetmax> <jetptmin> <checkjetid> <trkptmin> <gammaxi> <whichsys> <sysscalefactor> <inputfile>"
  exit 1
fi

# ./doff.sh 60 9999 30 0 1 0 0 0 DUMMYINPUT
# ./doff.sh DUMMY DUMMY DUMMY 0 1 0 0 0 DUMMYINPUT
runCmd="$HOME/code/scripts/myRun.sh"

set -x

inputPbPb="/mnt/hadoop/cms/store/user/tatar/GJT-out/skims/data_pbpb.root"
inputPP="/mnt/hadoop/cms/store/user/tatar/GJT-out/skims/data_pp.root"

centMins=(0  20  60 100);
centMaxs=(20 60 100 200);

phoPtMins=(60 80);
phoPtMaxs=(9999 9999);
jetPtMins=(30 40);

indicesCent=${!centMins[*]}
indicesPt=${!phoPtMins[*]}

for iCent in $indicesCent
do
  for iPt in $indicesPt
  do
    centMin=${centMins[iCent]}
    centMax=${centMaxs[iCent]}

    phoPtMin=${phoPtMins[iPt]}
    phoPtMax=${phoPtMaxs[iPt]}
    jetPtMin=${jetPtMins[iPt]}

    $runCmd ./ffgamma.exe $inputPbPb pbpbdata $centMin $centMax $phoPtMin $phoPtMax $jetPtMin recoreco $4 $5 $6 $7 $8 &> ffgamma_pbpbdata_${centMin}_${centMax}_${phoPtMin}_${phoPtMax}_${jetPtMin}_recoreco_${4}_${5}_${6}_${7}_${8}.log &
    wait

    $runCmd ./ffgamma.exe $inputPP ppdata $centMin $centMax $phoPtMin $phoPtMax $jetPtMin recoreco $4 $5 $6 $7 $8 &> ffgamma_ppdata_${centMin}_${centMax}_${phoPtMin}_${phoPtMax}_${jetPtMin}_recoreco_${4}_${5}_${6}_${7}_${8}.log &
    wait

    $runCmd hadd -f all_${phoPtMin}_${phoPtMax}_${jetPtMin}_gammaxi${6}.root pbpbdata_pbpbdata_*_*.root ppdata_ppdata_*_*.root &> hadd_all_${phoPtMin}_${phoPtMax}_${jetPtMin}_gammaxi${6}.log &
    wait

  done
done

rm pbpbdata_pbpbdata_*_*.root ppdata_ppdata_*_*.root 

