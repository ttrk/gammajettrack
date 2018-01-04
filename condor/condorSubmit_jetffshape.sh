#!/bin/bash

if [ $# -ne 14 ];
then
  echo "Usage: ./condor/condorSubmit_jetffshape.sh <skimFile> <sample> <hiBinMin> <hiBinMax> <phoetMin> <phoetMax> <jetptMin> <trkptMin> <gammaxi> <label> <systematics> <obs> <genlevel> <outputDir>"
  exit 1
fi

#g++ jetffshape.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o jetffshape.exe || exit 1
#echo "g++ jetffshape.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o jetffshape.exe || exit 1"

progFile="jetffshape.exe"

skimFile=$1
sample=$2
hiBinMin=$3
hiBinMax=$4
phoetMin=$5
phoetMax=$6
jetptMin=$7
trkptMin=$8
gammaxi=$9
label=${10}
systematics=${11}
obs=${12}
genlevel=${13}
outputDir=${14}

echo "skimFile : $skimFile"
echo "sample   : $sample"
echo "hiBinMin : $hiBinMin"
echo "hiBinMax : $hiBinMax"
echo "phoetMin : $phoetMin"
echo "phoetMax : $phoetMax"
echo "jetptMin : $jetptMin"
echo "trkptMin : $trkptMin"
echo "gammaxi  : $gammaxi"
echo "label    : $label"
echo "systematics : $systematics"
echo "obs : $obs"

outputName=$label"_"$sample"_"$genlevel"_"$phoetMin"_"$jetptMin"_"$gammaxi"_"$obs"_"$hiBinMin"_"$hiBinMax".root"

echo "outputDir  : $outputDir"
echo "outputName : $outputName"

# create the directories for condor submission and condor output files
baseDir="/work/"$USER"/gjt/"
timeNow=$(date +"%Y%m%d_%H%M%S")
outputNamePrefix="${outputName/.root/}"
submitDir=$baseDir"condorSubmissions/"$outputNamePrefix"/"$timeNow
condorLogsDir=$baseDir"condorLogs/"$outputNamePrefix"/"$timeNow
mkdir -p $submitDir
mkdir -p $condorLogsDir

echo "directory for condor submission : $submitDir" 
echo "directory for condor output     : $condorLogsDir" 

cp $progFile $submitDir
cp PbPb-weights.root $submitDir
cp pp-weights.root $submitDir
cp "condor/myRun.sh" $submitDir

## customizations for submit-hi2.mit.edu and submit.mit.edu machines ##
# proxy files start with "x509" and they are located under /tmp/ only.
proxyFilePath=$(find /tmp/ -maxdepth 1 -user $USER  -type f -name "x509*" -print | head -1)
proxyFile=$(basename $proxyFilePath)

srmPrefix="/mnt/hadoop/"
outputDirSRM=${outputDir#${srmPrefix}}

gfal-mkdir -p gsiftp://se01.cmsaf.mit.edu:2811/${outputDirSRM}
## customizations for submit-hi2.mit.edu and submit.mit.edu machines - END ##

# create the "submit description file"
cat > $submitDir/submit.condor <<EOF

Universe     = vanilla
Initialdir   = $submitDir
Notification = Error
Executable   = $submitDir/condorExecutable.sh
Arguments    = \$(Process) $progFile $skimFile $sample $hiBinMin $hiBinMax $phoetMin $phoetMax $jetptMin $trkptMin $gammaxi $label $systematics $obs $genlevel
GetEnv       = True
Output       = $condorLogsDir/\$(Process).out
Error        = $condorLogsDir/\$(Process).err
Log          = $condorLogsDir/\$(Process).log
Rank         = Mips
+AccountingGroup = "group_cmshi.$USER"
requirements = GLIDEIN_Site == "MIT_CampusFactory" && BOSCOGroup == "bosco_cmshi" && HAS_CVMFS_cms_cern_ch && BOSCOCluster == "ce03.cmsaf.mit.edu"
job_lease_duration = 240
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /tmp/$proxyFile,myRun.sh,$progFile,PbPb-weights.root,pp-weights.root

Queue $nJobs

EOF

# create the "executable" file
cat > $submitDir/condorExecutable.sh <<EOF
#!/bin/bash

progExe=\$2
skimFile=\$3
sample=\$4
hiBinMin=\$5
hiBinMax=\$6
phoetMin=\$7
phoetMax=\$8
jetptMin=\$9
trkptMin=\${10}
gammaxi=\${11}
label=\${12}
systematics=\${13}
obs=\${14}
genlevel=\${15}

outputDirTmp=$outputDir
outputTmp=$outputName

# setup grid proxy
export X509_USER_PROXY=\${PWD}/$proxyFile
# set hadoop directory path for grid tools
srmPrefix="/mnt/hadoop/"
outputDirSRM=\${outputDirTmp#\${srmPrefix}}

echo "###"
echo "host : \$(hostname)"
echo "PWD  : \$PWD"
echo "## voms-proxy-info --all ##"
voms-proxy-info --all
echo "###"
./myRun.sh "./"\$progExe \$skimFile \$sample \$hiBinMin \$hiBinMax \$phoetMin \$phoetMax \$jetptMin \$genlevel \$trkptMin \$gammaxi \$label \$systematics \$obs
echo "./myRun.sh "./"\$progExe \$skimFile \$sample \$hiBinMin \$hiBinMax \$phoetMin \$phoetMax \$jetptMin \$genlevel \$trkptMin \$gammaxi \$label \$systematics \$obs"

echo "## directory content ##"
ls -altrh
echo "## directory content - END ##"

set -x
mv -f \$outputTmp \$outputDirTmp

# $? is the exit status of last run command
if [ \$? -ne 0 ]; then
  gfal-copy -f file://\${PWD}/\${outputTmp} gsiftp://se01.cmsaf.mit.edu:2811/\${outputDirSRM}/\${outputTmp}

  if [ \$? -ne 0 ]; then
    srmcp -2 \$outputTmp gsiftp://se01.cmsaf.mit.edu:2811/\${outputDirSRM}/\${outputTmp}
  fi
fi
set +x

# delete all the files but the .out, .err, and .log files
ls | grep -v -e ".out" -e ".err" -e ".log" -e "_condor_stdout" -e "_condor_stderr" | xargs rm -rfv

echo "## directory content after deletion ##"
ls -altrh
echo "## directory content after deletion - END ##"

EOF

#cat $submitDir/submit.condor
echo "condor_submit $submitDir/submit.condor -name submit.mit.edu"
