#!/usr/local/bin/bash
#this is a general script that will take three bam files (kid, dad, mom) and generate a list of sites for which the kid is heterozygous
#this script calls ./check_snps_all.sh and ./trio_calls_vcfs/check_snps_all_trio.sh (requires a trio_config file in one-higher dir!) and get_base_counts.py but no changes are required to those scripts
#before using, change lines 6-11 and 41 (if exome data available)

trioname="ceu"
triolocation="/storage/CEUTrio/20110915/"
trio=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam')
trioshorthand=('878' '891' '892')   #kid first
chromosome="21"
#variable_site_file=ALL.chr${chromosome}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
variable_site_file="/storage/1kgenomes/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

pileupFile="chr${chromosome}_${trioshorthand[0]}.pileups"

exome=0     #if 1 make sure to input exome file in line 40
exome_file=CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.bam

#### Steup dir for two different callers
isOriginalCaller=false # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

if $isOriginalCaller;then
  workingDir="`pwd`/original/"
else
  workingDir="`pwd`/multiallelic/"
fi
if [ ! -d $workingDir ];then
  mkdir $workingDir
fi


####make a list to split up the chromosome for parallel runs
firstline=`samtools view ${triolocation}${trio[0]} ${chromosome} | head -1`
set -- $firstline
start=$4

#lastline=`samtools view ${triolocation}${trio[0]} ${chromosome} | tail -1`
#set -- $lastline
#end=$(($4+${#10}))

#Only use in example 
end=48119956

SplitFile="splitchr${chromosome}.split"

#num=$(($num1 + $num2))
#iterby=$((((${end}-${start}))/500))
#for i in `seq -f "%.0f" ${start} ${iterby} ${end}`;do echo "$i $(($i+((${iterby}-1)))) ${triolocation} ${trio[0]} ${trio[1]} ${trio[2]} ${trioname}";done > ${SplitFile}

#### Get trio and single calls
parallel_count=1
callsTrioDir="${workingDir}calls_trio_vcfs/" 
callsSingleDir="${workingDir}calls_single_vcfs/"

if [ ! -d ${callsTrioDir} ]; then
  mkdir ${callsTrioDir}
fi

if [ ! -d ${callsSingleDir} ]; then
  mkdir ${callsSingleDir}
fi

echo "=========="
#### In parallel, {} means input line. man parallel

####get trio calls - script for parallel
#cat ${SplitFile} | parallel -j $parallel_count ./calls_trio_vcfs/check_snps_all_trio.sh

#cat ${SplitFile}_test | xargs -I input ./check_snps_all_trio.sh input isOriginalCaller
cat ${SplitFile}_test |  parallel -j ${parallel_count} ./check_snps_all_trio.sh {} ${isOriginalCaller} ${callsTrioDir}

####get single calls - script for parallel
#cat ${SplitFile} | parallel -j $parallel_count ./check_snps_all.sh
cat ${SplitFile}_test | parallel -j ${parallel_count} ./check_snps_all.sh {} ${isOriginalCaller} ${callsSingleDir}

####get list of all vcf files
ls ${callsSingleDir}[1-9]*_single.vcf > ${workingDir}vcf_file_single_list.txt
ls ${callsTrioDir}[1-9]*_trio.vcf > ${workingDir}vcf_file_trio_list.txt



#### get all pileups. Reuse pileups for different bcftools mode
if [ ! -e ${pileupFile} ];then
    samtools mpileup -r ${chromosome} -f /storage/b37_reference/human_g1k_v37.fasta ${triolocation}${trio[0]} > ${pileupFile} 
fi


#### get all exome pileups. Reuse pileups for different bcftools mode
if [ $exome -eq 1 ];then
    if [ ! -e ${pileupExonFile} ];then
        samtools mpileup -r ${chromosome} -f /storage/b37_reference/human_g1k_v37.fasta ${triolocation}${exome_file} >  ${pileupExonFile}    
fi
fi

#sort hets and homos and use list of hets to get counts for those sites - for all trio and single calls, flags prev found as variable, exome
if [ $exome -eq 1 ];then
    python get_base_counts.py ${chromosome} ${trioshorthand[0]} ${variable_site_file} ${exome_file} ${workingDir}
else
    python get_base_counts.py ${chromosome} ${trioshorthand[0]} ${variable_site_file} 0 ${workingDir}
fi

#exit 10
#remove vcf files
#rm *single.vcf
#rm ./trio_calls_vcfs/*trio.vcf
