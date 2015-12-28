#!/usr/bin/env bash

#base_counts_hets_CEU10_878_byref_flag_filtered.txt       base_counts_homos_ref_CEU10_878_byref_flag_filtered.txt
#base_counts_homos_alt_CEU10_878_byref_flag_filtered.txt

# ln -sf /home/steven/Project_MDM/MiDiMu/R/run_gt_mdm.sh .

numTrial=1000

MdmFile="/home/steven/Project_MDM/MiDiMu/R/mdm.R"
GtMdm="/home/steven/Project_MDM/MiDiMu/R/search-gt-mdm.R"

#if [ -h $0 ]
#then
#	Orig=`readlink $0`
#	PWD=`dirname $Orig`
#else
#	Orig=$0
#	PWD=`pwd`
#fi
echo "Parameters: $0 $1 $#"

if [ $# -eq 1 ]  &&  [ "$1" == "D" ] || [ "$1" == "d" ] 
then
#	echo $1
	isDirty="d"
fi

#echo $Orig
#echo $PWD


#if [ ! -e "mdm.R" ] 
#then
#	ln -sf  $MdmFile .
#fi

#if [ ! -e "search-gt-mdm.R" ] 
#then
#	ln -sf $GtMdm .
#fi

Files=`ls | grep *hets*byref*`
echo "===== Estimate: $Files"
f=$Files
#for f in ${Files}
#do
	for numCom in $(seq 1 6)
	do
		echo "===== Process: $f ${numTrial} ${numCom}${isDirty}"
		#Rscript search-gt-mdm.R ${numTrial} ${numCom}${isDirty} $f
		Rscript ${GtMdm} ${numTrial} ${numCom}${isDirty} $f
	done
	

#done

