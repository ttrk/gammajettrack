#!/bin/bash

if [[ $# -ne 5 ]]; then
    echo "usage: ./skim-condor.sh [0/1/2/3] [input list] [output dir] [residuals] [mixing]"
    echo "[0: PbPb Data, 1: PbPb MC, 2: pp Data, 3: pp MC]"
    exit 1
fi

g++ ../photon_jet_track_skim.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o photon_jet_track_skim.exe
echo "g++ ../photon_jet_track_skim.C -Wall -Werror -O2 -Wno-narrowing `root-config --cflags --libs` -o photon_jet_track_skim.exe"

mkdir -p logs/
rm logs/*

mkdir -p $3

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
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = photon_jet_track_skim.exe,residuals.tgz,$2
#noop_job = !( stringListMember("\$(Process)","__FAILED__") )

Queue $JOBS
EOF

cat > skim.sh <<EOF
#!/bin/bash

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

echo ./photon_jet_track_skim.exe \$FILE \${1}.root \$JETALGO \$ISPP 1 \$MIXFILE
./photon_jet_track_skim.exe \$FILE \${1}.root \$JETALGO \$ISPP 1 \$MIXFILE

if [[ \$? -eq 0 ]]; then
    mv \${1}.root \${4}
fi

ls | grep -v .out | grep -v .err | grep -v .log | grep -v _condor_stdout | grep -v _condor_stderr | xargs rm -rf
EOF

condor_submit skim.condor
