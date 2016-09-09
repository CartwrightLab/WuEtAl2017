#!/usr/bin/env bash
# This is the setup script, variables and directories need to change accordingly
# This script calls pipelineCoreCHM1.sh, 
# which calls ./check_snps_all.sh and get_base_counts_chm1.py
#
# This will take one bam file,and generate a list of sites for which the kid is heterozygous


chm1Name="CHM1"
chromosome=$1  #21

ScriptHome="/home/steven/Project_MDM/MiDiMu/snps/"
parallel_count=20

if [ -z $chromosome ] ||  ([ $chromosome -ne 10 ] && [ $chromosome -ne 21 ]); then
    echo "ERROR: only chr 10 and 21 are supported now. Chromosome:${chromosome}"
    exit 2
fi

workingDir="/home/steven/Project_MDM/CHM1/${chm1Name}_C${chromosome}/"
if [ ! -d ${workingDir} ]; then
  mkdir -p ${workingDir}
fi


## sub chrome 10
if [ $chromosome -eq 10 ]; then
  chromosomePosition=":85534747-135534747"
#  echo "CP: ${chromosomePosition}"
fi


## VCF
variable_site_file="/storage/StevenStorage/CHM1/VCF/diploid_raw_UG_C${chromosome}.vcf"


chm1Location="/storage/chm1/bam/"
chm1Bam="SRR1283824.bam" 

echo "===== Setup ====="
echo "WorkingDir: $workingDir"
echo "VariableSiteFile: $variable_site_file"
echo "TrioLocation: ${chm1Location}"
#echo ${chm1Bam}

if [ -z ${chm1Location} ] || [ ! -d ${chm1Location} ];then
  echo "ERROR! Dir:${chm1Location} does not exits. chm1Name=${chm1Name}"
  exit 1
fi

isOriginalCaller=true #[true|false] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
exome=0    
exome_file=0 


# RUN

. ${ScriptHome}pipelineCoreCHM1.sh


