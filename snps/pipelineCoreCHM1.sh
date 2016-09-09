#!/usr/local/bin/bash
# This is the core script, should be called from another script (runCHM1_ChromY.sh) and variables are passed to this script.


# chm1Name="CHM1"
# chm1Location="/storage/CEUTrio/20130906/"
# variable_site_file="/storage/StevenStorage/CHM1/VCF/diploid_raw_UG_C${chromosome}.vcf"
# exome=0     #if 1 make sure to input exome file
# exome_file=0 
# parallel_count=20
# ScriptHome="/home/steven/Project_MDM/script/"



## FIX variable/scripts
run_check_snps_single="${ScriptHome}check_snps_all.sh"
python_get_base_counts="${ScriptHome}get_base_counts_chm1.py"

## END config, setup others (outfiles)

pileupFile="chr${chromosome}_${chm1Name}.pileups"

#### Steup dir for two different callers
##isOriginalCaller [true|false] # bcftools call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller

if $isOriginalCaller; then
    echo "===== bcftools original mode ====="
    workingDir="${workingDir}/original/"
else
    echo "===== bcftools multiallec mode ====="
    workingDir="${workingDir}/multiallelic/"
fi

if [ ! -d $workingDir ]; then
    mkdir $workingDir
fi

SplitFile="${workindDir}splitchr${chromosome}.split"

##### START #####
####make a list to split up the chromosome for parallel runs
firstline=`samtools view ${chm1Location}${chm1Bam} ${chromosome}${chromosomePosition} | head -1`
set -- $firstline
start=$4

lastline=`samtools view ${chm1Location}${chm1Bam} ${chromosome}${chromosomePosition} | tail -1`
set -- $lastline
end=$(($4+${#10}))


iterby=$((((${end}-${start}))/500))
for i in `seq -f "%.0f" ${start} ${iterby} ${end}`;do echo "$i $(($i+((${iterby}-1)))) ${chm1Location} ${chm1Bam} ${chm1Name}";done > ${SplitFile}



#### Get single calls
callsSingleDir="${workingDir}calls_single_vcfs/"


if [ ! -d ${callsSingleDir} ]; then
    mkdir ${callsSingleDir}
fi

#### print versions
echo "samtools --version"
samtools --version
echo "bcftools --version"
bcftools --version



#### In parallel, {} means input line. man parallel
####get single calls - script for parallel
echo "===== Check SNPS Single ====="

cat ${SplitFile} | parallel -j ${parallel_count} ${run_check_snps_single} ${chromosome} {} ${isOriginalCaller} ${callsSingleDir}


####get list of all vcf files
ls ${callsSingleDir}[1-9]*_single.vcf > ${workingDir}vcf_file_single_list.txt



#### get all pileups. Reuse pileups for different bcftools mode

if [[ ! ( -e ${pileupFile}  || -e ${pileupFile}.gz ) ]]; then
    echo "===== Create mpileup ======"
    samtools mpileup -r ${chromosome} -f /storage/reference_genomes/human/1k_genomes_phase1/human_g1k_v37.fasta ${chm1Location}${chm1Bam} > ${pileupFile} 
    gzip ${pileupFile}
else
    echo "FOUND ${pileupFile}.gz"
fi


### sort hets and homos and use list of hets to get counts for those sites - for all trio and single calls, flags prev found as variable, exome
if [ $exome -eq 1 ]; then
    python2 ${python_get_base_counts} ${chromosome} "${chm1Name}" ${variable_site_file} ${workingDir} ${pileupFile}.gz ${pileupExomeFile}
else
    echo "RUN: python2 ${python_get_base_counts} ${chromosome} "${chm1Name}" ${variable_site_file} ${workingDir} ${pileupFile}.gz 0"
    python2 ${python_get_base_counts} ${chromosome} "${chm1Name}" ${variable_site_file} ${workingDir} ${pileupFile}.gz 0
fi

# samtools mpileup -r 10 -f /storage/reference_genomes/human/1k_genomes_phase1/human_g1k_v37.fasta /storage/CEUTrio/20130906/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam > chr10_CEU13_878.pileups
# python2 get_base_counts.py 10 CEU13_878 /storage/StevenStorage/CEU/ALL.chr${chromosome}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf ~/Project_MDM/CEU/CEU13_C10/original/ chr10_CEU13_878.pileups 0


