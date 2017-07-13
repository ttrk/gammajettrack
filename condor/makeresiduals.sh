mkdir tmp
for i in `ls | grep -v .sh | grep -v .condor | grep -v logs | grep -v photon_jet_track_skim.exe | grep -v residuals | grep -v tmp `
do
    cp -r $i tmp
done
cd tmp
tar -czvf residuals.tgz *
mv residuals.tgz ..
cd ..
rm -rf tmp
