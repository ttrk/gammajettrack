pho=(40)
jet=(30)
ga=(0 1)
qg=(1 2)

for pt in 0
do  
    for g in `echo ${ga[*]}` 
    do
	for q in `echo ${qg[*]}` 
	do
	    hadd -f pyreweight_pbpb_${pho[$pt]}_1000_${jet[$pt]}_${q}_gammaxi${g}.root closure_pbpb_${pho[$pt]}_1000_${jet[$pt]}_${q}_gammaxi${g}.root v3syst/nominal/all_${pho[$pt]}_1000_${jet[$pt]}_gammaxi${g}.root  
	done
    done
done
