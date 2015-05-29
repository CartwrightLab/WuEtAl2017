#!/usr/local/bin/bash
# This is the core script, should be called from another script and variables are passed to this script.


# trioname="CEU13"
# triolocation="/storage/CEUTrio/20130906/"
# trio=('NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12892.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam')
# trioshorthand=('878' '891' '892')   #kid first
# chromosome="21"
# isOriginalCaller=false #[true|false] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

# variable_site_file=ALL.chr${chromosome}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
# variable_site_file="/storage/1kgenomes/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
# variable_site_file="/home/steven/Project_MDM/run/ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf"

# exome=0     #if 1 make sure to input exome file in line 40
# exome_file=0 

# parallel_count=3
#ScriptHome="/home/steven/Project_MDM/script/"



## FIX variable/scripts
run_check_snps_single="${ScriptHome}check_snps_all.sh"
run_check_snps_trio="${ScriptHome}check_snps_all_trio.sh"
python_get_base_counts="${ScriptHome}get_base_counts.py"

## END config, setup others (outfiles)

pileupFile="chr${chromosome}_${trioname}_${trioshorthand[0]}.pileups"
pileupExomeFile="chr${chromosome}Ex_${trioname}_${trioshorthand[0]}.pileups"


#### Steup dir for two different callers
##isOriginalCaller [true|false] # bcftools call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

if $isOriginalCaller;then
  echo "===== bcftools original mode ====="
  workingDir="${workingDir}/original/"
else
  echo "===== bcftools multiallec mode ====="
  workingDir="${workingDir}/multiallelic/"
fi
if [ ! -d $workingDir ];then
  mkdir $workingDir
fi

SplitFile="${workindDir}splitchr${chromosome}.split"

####make a list to split up the chromosome for parallel runs
firstline=`samtools view ${triolocation}${trio[0]} ${chromosome}${chromosomePosition} | head -1`
set -- $firstline
start=$4

lastline=`samtools view ${triolocation}${trio[0]} ${chromosome}${chromosomePosition} | tail -1`
set -- $lastline
end=$(($4+${#10}))

#end=48119956 ##ONLY use this is the example to speed things up
#echo $firstline
#echo $lastline
#echo $start
#echo $end
#exit 5

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

cat ${SplitFile} |  parallel -j ${parallel_count} ${run_check_snps_trio} ${chromosome} {} ${isOriginalCaller} ${callsTrioDir}
#cat ${SplitFile}_test | xargs -I input ./check_snps_all_trio.sh input isOriginalCaller


####get single calls - script for parallel
echo "===== Check SNPS Single ====="

cat ${SplitFile} | parallel -j ${parallel_count} ${run_check_snps_single} ${chromosome} {} ${isOriginalCaller} ${callsSingleDir}
#cat ${SplitFile} | parallel -j ${parallel_count} ./check_snps_all.sh {} ${isOriginalCaller} ${callsSingleDir}

####get list of all vcf files
ls ${callsSingleDir}[1-9]*_single.vcf > ${workingDir}vcf_file_single_list.txt
ls ${callsTrioDir}[1-9]*_trio.vcf > ${workingDir}vcf_file_trio_list.txt



#### get all pileups. Reuse pileups for different bcftools mode
if [ ! -e ${pileupFile} ];then
    echo "===== Create mpileup ======"
    samtools mpileup -r ${chromosome} -f /storage/reference_genomes/human/1k_genomes_phase1/human_g1k_v37.fasta ${triolocation}${trio[0]} > ${pileupFile} 
fi


#### get all exome pileups. Reuse pileups for different bcftools mode
if [ $exome -eq 1 ];then
    if [ ! -e ${pileupExomeFile} ];then
        echo "===== Create mpileup exome ======"
        samtools mpileup -r ${chromosome} -f /storage/reference_genomes/human/1k_genomes_phase1/human_g1k_v37.fasta ${triolocation}${exome_file} >  ${pileupExomeFile}    
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
