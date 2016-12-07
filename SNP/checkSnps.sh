#!/usr/bin/bash
# Call SNP for trio or a single individual using Samtools and BCFtools
### Requirement
# - [Samtools](http://www.htslib.org/) version 1.2 or later
# - [BCFtools](http://www.htslib.org/) version 1.2 or later
#
# Trio: 9 parameters required.
# Usage: check_snps.sh chromosome isOriginalCaller outfileDir refGenome
#   infilePrefix infileChild infileParent1 infileParent2
#   "startPos endPos"
#
# Single: 6 parameters required.
# Usage: check_snps.sh chromosome isOriginalCaller outfileDir refGenome infile
#   "startPos endPos"

DEBUG=false
parseStartEnd (){
  if [ $# -eq 2 ]; then
    startpos=$1
    endpos=$2
    if [ $endpos -lt $startpos ];then
      echo -e "Error: starting position is greater than ending position ($startpos > $endpos).\nExit!"
      exit 1
    fi

  else
    echo -e "Error: Cannot pares position information. ($@)\nExit!"
    echo "Syntax: \"startPosition endPosition\""
    exit 1
  fi
}

chromosome=$1
isOriginalCaller=$2 #[1|0]
outfileDir=$3
refGenome=$4

if [ $# -eq 9 ]; then
  isTrio=1
  infile=($5$6 $5$7 $5$8) #child parent1 parent2
  parseStartEnd $9

elif [[ $# -eq 6 ]]; then
  isTrio=0
  infile=$5
  parseStartEnd $6

else
  echo -e "Error: incorrect numebr of parameters ($#). Exit!"
  echo -e "\tTrio: 9 parameters required.\n\tUsage: checkSnps.sh chromosome isOriginalCaller outfileDir refGenome\n\t\tinfilePrefix infileChild infileParent1 infileParent2\n\t\t\"startPos endPos\"\n"
  echo -e "\tSingle: 6 parameters required.\n\tUsage: checkSnps.sh chromosome isOriginalCaller outfileDir refGenome infile\n\t\t\"startPos endPos\""
  echo "Command: $@"
  exit 1
fi

if [ $isOriginalCaller -eq 1 ]; then
  mode="-c" #-c, --consensus-caller	bcftools
else
  mode="-m" #-m, --multiallelic-caller	bcftools
fi


outfile="${outfileDir}${startpos}.vcf"
if [ -e $outfile ]; then
  lastline=`tail -1 ${outfile}`  #check for restarting - if so get last line of vcf file
  set -- $lastline
  if [ "$2" != "POS" ] &&  [ $2 -gt $startpos ] && [ $2 -lt $endpos ];then    #get last position tested
  #    echo "reset $startpos to $2+1 END:$endpos"
      startpos=$(($2+1)) #start with +1 position, got result from this position already
      outfile="${outfileDir}${startpos}.vcf"
  fi
fi

if [ ! -e $outfile ];then
  if [ $isTrio -eq 1 ]; then
    echo "===== Check SNPS Trio: starting ${startpos} ====="
    ( samtools mpileup -t SP -u -r ${chromosome}:${startpos}-${endpos} -f ${refGenome} ${infile[0]} ${infile[1]} ${infile[2]} | bcftools call ${mode} - > ${outfile} )
  else
    echo "===== Check SNPS Single: starting ${startpos} ====="
    ( samtools mpileup -t SP -u -r ${chromosome}:${startpos}-${endpos} -f   ${refGenome} ${infile[0]} | bcftools call ${mode} - > ${outfile} )
  fi
fi


exit 0
