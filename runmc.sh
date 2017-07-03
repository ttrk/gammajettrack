pho=(40) # pho et bins
jet=(30) # jet pt bins
ga=(0)      # 0 = xijet , 1 = xigamma
qg=(0 1 2)  # 0 = inclusive , 1 = quark only, 2 = gluon only

for pt in 0 
do  
    for g in `echo ${ga[*]}` 
    do
	for q in `echo ${qg[*]}` 
	do
	    if [[ $1 -eq 999 ]] # run jet loop only
	    then
		./domcjetpt.sh ${pho[$pt]} 1000 ${jet[$pt]} ${q} 1 ${g} recoreco gengen0 
	    else              # default
		./doclosure.sh ${pho[$pt]} 1000 ${jet[$pt]} ${q} 1 ${g} recoreco gengen0 
	    fi
	done
    done
done
