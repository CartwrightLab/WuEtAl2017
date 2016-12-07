#!/usr/bin/bash
### Requirement
# - [Samtools](http://www.htslib.org/) version 1.2 or later
# - [BCFtools](http://www.htslib.org/) version 1.2 or later
# - [Python 2.7](https://www.python.org/)
# - **Recommendation:** [GNU Parallel](https://www.gnu.org/software/parallel/)


## Requires the following variables to setup
datasetName="CEU13_example" #Can be any name
chromosome=21 #Which chromosome would you like to perforem the analyses on
chromosomePosition=":9411180-9421180"  #Perfrom analyses on a specific region e.g. chromosomePosition=":10000-20000"
##chromosomePosition="" #Blank: perform analyses on the whole chromosame.

workingDir="./data/${datasetName}_C${chromosome}/"  #Location of the output files
scriptDir="./"  #Location of the scripts

variableSiteFile="./data/snp_chr21_v3_partial.vcf.gz" #location of known snp files
refGenome="./data/human_g1k_v37_c21_partial.fasta" #Location of the reference file
dataDir="./data/"   #Location of the datafile
dataFiles=('CEU13_NA12878.bam' 'CEU13_NA12892.bam' 'CEU13_NA12891.bam') #File names of each individual

parallelCount=3 #Number of parallel jobs
numSplitFiles=5 #Number of split files
isTrio=1 #1:trio 0:single caller
isOriginalCaller=1 #[1|0] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
exome=0 #1: turn on exome analyses
exome_file=0 #Location of the exome files

### End of user parameters inputs



. ${scriptDir}setupParameters.sh
## Done with setup
if [ ! -f ${splitFile} ]; then
  echo "===== Build split files for parallel runs ====="
  firstline=`samtools view ${dataDir}${dataFiles[0]} ${chromosome}${chromosomePosition} | head -1`
  set -- $firstline
  start=$4

  lastline=`samtools view ${dataDir}${dataFiles[0]} ${chromosome}${chromosomePosition} | tail -1`
  set -- $lastline
  end=$(($4+${#10}))

  iterby=$((((${end}-${start}))/${numSplitFiles}))
  for i in `seq -f "%.0f" ${start} ${iterby} ${end}`; do echo "$i $(($i+((${iterby}-1))))"; done > ${splitFile}
else
  echo "===== FOUND ${splitFile} and reuse this. ====="
fi


## Run Check Snps Script
echo "===== Start SNPS Single ====="
cat ${splitFile} |  ${runCommand} ${runCheckSnpsScript} ${chromosome} ${isOriginalCaller} ${callsSingleDir} ${refGenome} ${dataDir}${dataFiles[0]} {}

ls ${callsSingleDir}[1-9]*.vcf > ${workingDir}vcf_file_single_list.txt

if [ $isTrio == 1 ];then
  echo "===== Start SNPS Trio ====="
  cat ${splitFile} |  ${runCommand} ${runCheckSnpsScript} ${chromosome} ${isOriginalCaller} ${callsTrioDir} ${refGenome} ${dataDir} ${dataFiles[0]} ${dataFiles[1]} ${dataFiles[2]} {}

  ls ${callsTrioDir}[1-9]*.vcf > ${workingDir}vcf_file_trio_list.txt
fi

## Create or reuse pileup
if [ !  -f ${pileupFile}.gz  ]; then
    echo "===== Create mpileup ======"
    samtools mpileup -r ${chromosome} -f ${refGenome} ${dataDir}${dataFiles[0]} > ${pileupFile}
    gzip ${pileupFile}
else
    echo "===== FOUND ${pileupFile}.gz and reuse this. ====="
fi

if [ $exome -eq 1 ]; then
    if [[ ! ( -e ${pileupExomeFile}  || -e ${pileupExomeFile}.gz ) ]]; then
        echo "===== Create mpileup exome ======"
        samtools mpileup -r ${chromosome} -f ${refGenome} ${dataDir}${exome_file} >  ${pileupExomeFile}
        gzip ${pileupExomeFile}
    else
        echo "===== FOUND ${pileupExomeFile}.gz and reuse this. ====="
    fi
else
  pileupExomeFile=0
fi

## sort hets and homos and use list of hets to get counts for those sites - for all trio and single calls, flags prev found as variable, exome
python2 ${pythonGetBaseCounts} ${chromosome} ${datasetName} ${variableSiteFile} ${workingDir} ${pileupFile}.gz ${pileupExomeFile}
