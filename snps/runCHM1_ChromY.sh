#!/usr/bin/env bash
#this is the setup script, call pipelineCoreCHM1.sh
# that will take one bam file,and generate a list of sites for which the kid is heterozygous
#this script calls ./check_snps_all.sh and get_base_counts.py

chm1Name="CHM1" #"CHM1" FIXED
chromosome=$1  #21

ScriptHome="/home/steven/Project_MDM/MiDiMu/snps/"
parallel_count=20

#if [ ! -z $chromosome ];then
#    echo "No chromosome seleted $chromosome"
#fi

if [ -z $chromosome ] ||  ([ $chromosome -ne 10 ] && [ $chromosome -ne 21 ]) ; then
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
variable_site_file="/storage/StevenStorage/CEU/ALL.chr${chromosome}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
variable_site_file="/storage/StevenStorage/CHM1/VCF/diploid_raw_UG_C${chromosome}.vcf"


## FIX setting
#trioshorthand=('878' '891' '892')   #kid first

chm1Location="/storage/chm1/bam/"
chm1Bam="SRR1283824.bam" #('NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12892.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam')

echo "===== Setup ====="
echo "WorkingDir: $workingDir"
echo "VariableSiteFile: $variable_site_file"
echo "TrioLocation: ${chm1Location}"
#echo ${chm1Bam}

if [ -z ${chm1Location} ] || [ ! -d ${chm1Location} ]
then
  echo "ERROR! Dir:${chm1Location} does not exits. chm1Name=${chm1Name}"
  exit 1
fi

isOriginalCaller=true #[true|false] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
exome=0    
exome_file=0 


# RUN

. ${ScriptHome}pipelineCoreCHM1.sh


