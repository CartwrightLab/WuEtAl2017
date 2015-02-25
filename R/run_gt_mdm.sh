#!/bin/sh

#base_counts_hets_CEU10_878_byref_flag_filtered.txt       base_counts_homos_ref_CEU10_878_byref_flag_filtered.txt
#base_counts_homos_alt_CEU10_878_byref_flag_filtered.txt
#PWD=`pwd`
#MdmFile="/home/steven/Postdoc2/Project_MDM/RCode/mdm.R"
#GtMdm="/home/steven/Postdoc2/Project_MDM/RCode/search-gt-mdm.R"

if [ -h $0 ]
then
	Orig=`readlink $0`
	PWD=`dirname $Orig`
else
	Orig=$0
	PWD=`pwd`
fi

if [ $1 -eq "D" | $1 -eq "d" ]
then
	isDirty="d"
fi

MdmFile="${PWD}/mdm.R"
GtMdm="${PWD}/search-gt-mdm.R"
echo $Orig
echo $PWD

numTrial=20

if [ ! -e "mdm.R" ] 
then
	ln -sf  $MdmFile .
fi

if [ ! -e "search-gt-mdm.R" ] 
then
	ln -sf $GtMdm .
fi

Files=`ls | grep byref`
for f in ${Files}
do
	for numCom in $(seq 1 6)
	do
#		echo $f $numCom
		Rscript search-gt-mdm.R ${numTrial} ${numCom}${isDirty} $f
	done
	

done

