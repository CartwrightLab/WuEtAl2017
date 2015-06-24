#!/bin/sh
chromosome=$1
line=$2     #range in chromosome
isOriginalCaller=$3
outfileDir=$4

set -- $line
startpos=$1
endpos=$2
infile=($3$4 $3$5 $3$6) #child parent1 parent2

if $isOriginalCaller; then
  mode="-c" #-c, --consensus-caller	bcftools
else
  mode="-m" #-m, --multiallelic-caller	bcftools
fi

outfile="${outfileDir}${startpos}_trio.vcf"
[ -e $outfile ] && lastline=`tail -1 ${outfile}` || lastline='1 1'  #check for restarting - if so get last line of vcf file otherwise use a pretend line
set -- $lastline        
if [[ $2 -gt $startpos ]];then    #get last position tested
#    echo "reset $startpos to $2+1 END:$endpos"
    startpos=$(($2+1)) #start with +1 position, got result from this position already
    outfile="${outfileDir}${startpos}_trio.vcf"
#else
#echo "Skip T1"
fi


if [ $endpos -gt $startpos ];then
  ( samtools mpileup -t SP -u -r ${chromosome}:${startpos}-${endpos} -f /storage/reference_genomes/human/1k_genomes_phase1/human_g1k_v37.fasta ${infile[0]} ${infile[1]} ${infile[2]} | bcftools call ${mode} -C trio -S trio_config_CEU.txt -> ${outfile} )
#else
#echo "Skip T2"
fi

#Alternative use inline parameters instead of config file
#-s NA19240 NA19238 NA19239

# NEW Version
# bcftools --version
## bcftools 1.1
## Using htslib 1.1
## Copyright (C) 2014 Genome Research Ltd.
# samtools --version
## samtools 1.1
## Using htslib 1.1
## Copyright (C) 2014 Genome Research Ltd.

#( samtools mpileup -t SP -u -r 21:${startpos}-${endpos} -f /storage/b37_reference/human_g1k_v37.fasta /storage/yri_trio/NA19240.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19238.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19239.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam | bcftools call -c -C trio -S trio_config_yri.txt -> ${startpos}.vcf )
#note starting new file

#-c, --consensus-caller         the original samtools/bcftools calling method (conflicts with -m)
# -m, --multiallelic-caller     alternative modelfor multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c)



#Old version samtools/bcftools
# samtools
# -D Deprecated => -t SP
# bcftools view => bcftools call i

# bcftools old parameters
#-c 	Call variants using Bayesian inference. This option automatically invokes option -e.
#-g	Call per-sample genotypes at variant sites (force -c)
#-T STR	Enable pair/trio calling. For trio calling, option -s is usually needed to be applied to configure the trio members and their ordering. In the file supplied to the option -s, the first sample must be the child, the second the father and the third the mother. The valid values of STR are ‘pair’, ‘trioauto’, ‘trioxd’ and ‘trioxs’, where ‘pair’ calls differences between two input samples, and ‘trioxd’ (‘trioxs’) specifies that the input is from the X chromosome non-PAR regions and the child is a female (male). [null]
#-s FILE	List of samples to use. The first column in the input gives the sample names and the second gives the ploidy, which can only be 1 or 2. When the 2nd column is absent, the sample ploidy is assumed to be 2. In the output, the ordering of samples will be identical to the one in FILE. [null]

#( samtools mpileup -uD -r 21:${startpos}-${endpos} -f /storage/b37_reference/human_g1k_v37.fasta /storage/yri_trio/NA19240.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19238.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam /storage/yri_trio/NA19239.chrom21.ILLUMINA.bwa.YRI.high_coverage.20100311.bam | bcftools view -cgT trioauto -s ../../trio_config_yri.txt -> ${startpos}.vcf )      #note starting new file


exit 0

#cat ../splitchr21.txt | parallel -j 60 ./check_snps_all_trio.sh
