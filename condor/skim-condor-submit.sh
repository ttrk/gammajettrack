#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "usage: ./skim-condor-submit.sh [0/1/2/3] [input list] [output dir] [residuals] [mixing]"
    echo "[0: PbPb Data, 1: PbPb MC, 2: pp Data, 3: pp MC]"
    exit 1
fi

g++ ../photon_jet_track_skim.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o photon_jet_track_skim.exe
echo "g++ ../photon_jet_track_skim.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o photon_jet_track_skim.exe"

PROXYFILE=$(ls /tmp/ -lt | grep $USER | grep -m 1 x509 | awk '{print $NF}')

SRM_PREFIX="/mnt/hadoop/"
SRM_PATH=${3#${SRM_PREFIX}}

gfal-mkdir -p gsiftp://se01.cmsaf.mit.edu:2811/${SRM_PATH}

[ -d "logs" ] && rm -r logs
mkdir -p logs/

JOBS=$(cat $2 | wc -l)

cat > skim.condor <<EOF
Universe     = vanilla
Initialdir   = $PWD/
Notification = Error
Executable   = $PWD/skim.sh
Arguments    = \$(Process) $1 $2 $3 $4 $5
GetEnv       = True
Output       = $PWD/logs/\$(Process).out
Error        = $PWD/logs/\$(Process).err
Log          = $PWD/logs/\$(Process).log
Rank         = Mips
+AccountingGroup = "group_cmshi.$(whoami)"
requirements = GLIDEIN_Site == "MIT_CampusFactory" && BOSCOGroup == "bosco_cmshi" && HAS_CVMFS_cms_cern_ch && BOSCOCluster == "ce03.cmsaf.mit.edu"
job_lease_duration = 240
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /tmp/$PROXYFILE,photon_jet_track_skim.exe,$4,$2
#noop_job = !( stringListMember("\$(Process)","__FAILED__") )

Queue $JOBS
EOF

cat > skim.sh <<EOF
#!/bin/bash

ls
echo \$(whoami)
echo \$HOSTNAME

# setup grid proxy
export X509_USER_PROXY=\${PWD}/$PROXYFILE

# set hadoop directory path for grid tools
SRM_PREFIX="/mnt/hadoop/"
SRM_PATH=\${4#\${SRM_PREFIX}}

# extract corrections and other files needed to run skim
tar -xzvf \$5

FILE=\$(head -n\$((\$1+1)) \$3 | tail -n1)
case \$2 in
    0)
        JETALGO=akPu3PFJetAnalyzer
        ISPP=0
        NMIX=\$((\$1%4))
        NMIX=\$((\$NMIX+1))
        MIXFILE=\$(head -n\${NMIX} \$6 | tail -n1)
        ;;
    1)
        JETALGO=akPu3PFJetAnalyzer
        ISPP=0
        MIXFILE=/mnt/hadoop/cms/store/user/biran/photon-jet-track/PbPb-MB-Hydjet-Cymbal-170331.root
        ;;
    2)
        JETALGO=ak3PFJetAnalyzer
        ISPP=1
        MXIFILE=""
        ;;
    3)
        JETALGO=ak3PFJetAnalyzer
        ISPP=1
        MIXFILE=""
        ;;
    *)
        echo "bad option"
        exit 1
        ;;
esac

set -x

echo ./photon_jet_track_skim.exe \$FILE \${1}.root \$JETALGO \$ISPP 1 \$MIXFILE
./photon_jet_track_skim.exe \$FILE \${1}.root \$JETALGO \$ISPP 1 \$MIXFILE

if [[ \$? -eq 0 ]]; then
    mv \${1}.root \${4}

    if [[ \$? -ne 0 ]]; then
        gfal-copy file://\$PWD/\${1}.root gsiftp://se01.cmsaf.mit.edu:2811/\${SRM_PATH}/\${1}.root

        if [[ \$? -ne 0 ]]; then
            srmcp -2 \${1}.root gsiftp://se01.cmsaf.mit.edu:2811/\${SRM_PATH}/\${1}.root
        fi
    fi
fi

ls | grep -v .out | grep -v .err | grep -v .log | grep -v _condor_stdout | grep -v _condor_stderr | xargs rm -rf
EOF

condor_submit skim.condor -pool submit.mit.edu:9615 -name submit.mit.edu -spool
