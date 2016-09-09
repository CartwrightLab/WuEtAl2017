#!/usr/bin/env bash
## Run EM on files with the following pattern *hets*byref*
## e.g. base_counts_hets_CHM1_byref_flag_filtered.txt

numTrial=1000

MdmFile="./mdm.R"
GtMdm="./search_gt_chm1_mdm.R"

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

