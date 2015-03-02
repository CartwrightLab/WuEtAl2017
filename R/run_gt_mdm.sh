#!/bin/sh

#base_counts_hets_CEU10_878_byref_flag_filtered.txt       base_counts_homos_ref_CEU10_878_byref_flag_filtered.txt
#base_counts_homos_alt_CEU10_878_byref_flag_filtered.txt


if [ -h $0 ]
then
	Orig=`readlink $0`
	PWD=`dirname $Orig`
else
	Orig=$0
	PWD=`pwd`
fi
echo $1 $#
if [[ $# -eq 1  &&  $1 == "D" || $1 == "d" ]]
then
	echo $1
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

Files=`ls | grep *hets*byref*`
for f in ${Files}
do
	for numCom in $(seq 1 6)
	do
		echo "$f ${numTrial} ${numCom}${isDirty}"
		Rscript search-gt-mdm.R ${numTrial} ${numCom}${isDirty} $f
	done
	

done

