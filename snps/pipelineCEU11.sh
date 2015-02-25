#!/usr/local/bin/bash
#this is a general script that will take three bam files (kid, dad, mom) and generate a list of sites for which the kid is heterozygous
#this script calls ./check_snps_all.sh and ./check_snps_all_trio.sh and get_base_counts.py but no changes are required to those scripts
#before using, change lines 6-21 

trioname="CEU11"
triolocation="/storage/CEUTrio/20110915/"
trio=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam')
trioshorthand=('878' '891' '892')   #kid first
chromosome="21"
isOriginalCaller=true #[true|false] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

#variable_site_file=ALL.chr${chromosome}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#variable_site_file="/storage/1kgenomes/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
variable_site_file="/home/steven/Project_MDM/run/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf"

exome=1     #if 1 make sure to input exome file in line 40
exome_file=CEUTrio.HiSeq.WEx.b37_decoy.NA12878.clean.dedup.recal.bam

ScriptHome="/home/steven/Project_MDM/script/snps/"
run_check_snps_single="${ScriptHome}check_snps_all.sh"
run_check_snps_trio="${ScriptHome}check_snps_all_trio.sh"
python_get_base_counts="${ScriptHome}get_base_counts.py"

## END config, setup others

pileupFile="chr${chromosome}_${trioname}_${trioshorthand[0]}.pileups"
pileupExomeFile="chr${chromosome}Ex_${trioname}_${trioshorthand[0]}.pileups"

parallel_count=5
#### Steup dir for two different callers
##isOriginalCaller [true|false] # bcftools call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

if $isOriginalCaller;then
  echo "===== bcftools original mode ====="
  workingDir="`pwd`/original/"
else
  echo "===== bcftools multiallec mode ====="
  workingDir="`pwd`/multiallelic/"
fi
if [ ! -d $workingDir ];then
  mkdir $workingDir
fi

SplitFile="${workindDir}splitchr${chromosome}.split"

####make a list to split up the chromosome for parallel runs
firstline=`samtools view ${triolocation}${trio[0]} ${chromosome} | head -1`
set -- $firstline
start=$4

lastline=`samtools view ${triolocation}${trio[0]} ${chromosome} | tail -1`
set -- $lastline
end=$(($4+${#10}))

#end=48119956 ##ONLY use this is the example to speed things up


##num=$(($num1 + $num2))

iterby=$((((${end}-${start}))/500))
for i in `seq -f "%.0f" ${start} ${iterby} ${end}`;do echo "$i $(($i+((${iterby}-1)))) ${triolocation} ${trio[0]} ${trio[1]} ${trio[2]} ${trioname}";done > ${SplitFile}



#### Get trio and single calls
callsTrioDir="${workingDir}calls_trio_vcfs/" 
callsSingleDir="${workingDir}calls_single_vcfs/"

if [ ! -d ${callsTrioDir} ]; then
  mkdir ${callsTrioDir}
fi

if [ ! -d ${callsSingleDir} ]; then
  mkdir ${callsSingleDir}
fi

#### In parallel, {} means input line. man parallel

####get trio calls - script for parallel
echo "===== Check SNPS Trio ====="

cat ${SplitFile} |  parallel -j ${parallel_count} ${run_check_snps_trio} {} ${isOriginalCaller} ${callsTrioDir}
#cat ${SplitFile}_test | xargs -I input ./check_snps_all_trio.sh input isOriginalCaller


####get single calls - script for parallel
echo "===== Check SNPS Single ====="

cat ${SplitFile} | parallel -j ${parallel_count} ${run_check_snps_single} {} ${isOriginalCaller} ${callsSingleDir}
#cat ${SplitFile} | parallel -j ${parallel_count} ./check_snps_all.sh {} ${isOriginalCaller} ${callsSingleDir}

####get list of all vcf files
ls ${callsSingleDir}[1-9]*_single.vcf > ${workingDir}vcf_file_single_list.txt
ls ${callsTrioDir}[1-9]*_trio.vcf > ${workingDir}vcf_file_trio_list.txt



#### get all pileups. Reuse pileups for different bcftools mode
if [ ! -e ${pileupFile} ];then
    echo "===== Create mpileup ======"
    samtools mpileup -r ${chromosome} -f /storage/b37_reference/human_g1k_v37.fasta ${triolocation}${trio[0]} > ${pileupFile} 
fi


#### get all exome pileups. Reuse pileups for different bcftools mode
if [ $exome -eq 1 ];then
    if [ ! -e ${pileupExomeFile} ];then
        echo "===== Create mpileup exome ======"
        samtools mpileup -r ${chromosome} -f /storage/b37_reference/human_g1k_v37.fasta ${triolocation}${exome_file} >  ${pileupExomeFile}    
fi
fi

### sort hets and homos and use list of hets to get counts for those sites - for all trio and single calls, flags prev found as variable, exome
if [ $exome -eq 1 ];then
    python ${python_get_base_counts} ${chromosome} "${trioname}_${trioshorthand[0]}" ${variable_site_file} ${workingDir} ${pileupFile} ${pileupExomeFile}
else
    python ${python_get_base_counts} ${chromosome} "${trioname}_${trioshorthand[0]}" ${variable_site_file} ${workingDir} ${pileupFile} 0
fi


#remove vcf files
#rm *single.vcf
#rm ./trio_calls_vcfs/*trio.vcf
