#/bin/sh

line=$1     #range in chromosome
isOriginalCaller=$2
outfileDir=$3

set -- $line
startpos=$1
endpos=$2
infile=($3$4 $3$5 $3$6) #child parent1 parent2

if $isOriginalCaller; then
  mode="-c" #-c, --consensus-caller	bcftools
else
  mode="-m" #-m, --multiallelic-caller	bcftools
fi

outfile="${outfileDir}${startpos}_single.vcf"
[ -e $outfile ] && lastline=`tail -1 ${outfile}` || lastline='1 1'  #check for restarting - if so get last line of vcf file otherwise use a pretend line
set -- $lastline        
if [ "$2" -gt "$startpos" ];then    #get last position tested
#    echo "reset $startpos to $2+1 END:$endpos"
    startpos=$(($2+1)) #start with +1 position, got result from this position already
    outfile="${outfileDir}${startpos}_single.vcf"
fi


if [ $endpos -gt $startpos ];then
  ( samtools mpileup -t SP -u -r 21:${startpos}-${endpos} -f /storage/b37_reference/human_g1k_v37.fasta ${infile[0]} | bcftools call -c -> ${outfile} )      #note starting new file
fi

##orig
#( samtools mpileup -uD -r 21:${startpos}-${endpos} -f /storage/b37_reference/human_g1k_v37.fasta /storage/CEUTrio/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam | bcftools view -cg -> ${startpos}_878.vcf )      #note starting new file

#  samtools mpileup -t SP -u -r 21:${startpos}-${endpos} -f /storage/b37_reference/human_g1k_v37.fasta /storage/yri_trio/NA19240.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19238.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19239.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam | bcftools call -c -C trio -S trio_config_yri.txt -> ${startpos}.vcf )



# NEW Version
# bcftools --version
## bcftools 1.1
## Using htslib 1.1
## Copyright (C) 2014 Genome Research Ltd.
# samtools --version
## samtools 1.1
## Using htslib 1.1
## Copyright (C) 2014 Genome Research Ltd.

#bcftolos
#-c, --consensus-caller		the original samtools/bcftools calling method (conflicts with -m)
# -m, --multiallelic-caller	alternative modelfor multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c)

#cat splitchr21.txt | parallel -j 60 ./check_snps_all.sh
