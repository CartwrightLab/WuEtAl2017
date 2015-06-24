#!/usr/bin/env bash
#this is the setup script, call pipelineCore.sh
# that will take three bam files (kid, dad, mom) and generate a list of sites for which the kid is heterozygous
#this script calls ./check_snps_all.sh and ./check_snps_all_trio.sh and get_base_counts.py

trioname="CEU${1}" #"CEU13"
chromosome=$2 #21

ScriptHome="/home/steven/Project_MDM/MiDiMu/snps/"
parallel_count=20

workingDir="/home/steven/Project_MDM/CEU/${trioname}_C${chromosome}/"
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


## FIX setting
trioshorthand=('878' '891' '892')   #kid first

if [ ${trioname} == "CEU10" ] 
then
  triolocation="/storage/CEUTrio/20100311/"
  trio=('NA12878.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100311.bam' 'NA12892.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam' 'NA12891.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam')
  if [ $chromosome -eq 10 ]; then
    trio=('NA12878.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100311.bam' 'NA12891.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam' 'NA12892.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam')
  fi

elif [ ${trioname} == "CEU11" ]
then
  triolocation="/storage/CEUTrio/20110915/"
  trio=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam')
elif [ ${trioname} == "CEU12" ]
then
  triolocation="/storage/CEUTrio/20120117/"
  trio=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.20120117.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.20120117.bam')
elif [ ${trioname} == "CEU13" ]
then
  triolocation="/storage/CEUTrio/20130906/"
  trio=('NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12892.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam')
fi

echo "===== Setup ====="
echo "WorkingDir: $workingDir"
echo "VariableSiteFile: $variable_site_file"
echo "TrioLocation: ${triolocation}"
#echo ${trio[0]} ${trio[1]} ${trio[2]}

if [ -z ${triolocation} ] || [ ! -d ${triolocation} ]
then
  echo "ERROR! Dir:${triolocation} does not exits. trioname=${trioname}"
  exit 1
fi

isOriginalCaller=true #[true|false] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
exome=0    
exome_file=0 


# RUN

. ${ScriptHome}pipelineCore.sh


