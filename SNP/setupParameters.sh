#!/usr/bin/env bash
## Check program and parameters
## Check the availability of samtools, brftools, and parallel
## If parallel is installed, run it with ${parallelCount} jobs
## Check if ${dataDir} and ${scriptDir} exist
## Setup variables for the pipeline

## This script setup the following variables
## $runCommand           = command used to run check_snps.sh.
## $workingDir           = working directory for storing final results.
## $callsSingleDir       = directory to store vcf fils for single calls.
## callsTrioDir          = directory to store vcf fils for trio calls. Only set if isTrue==1.
## $pythonGetBaseCounts  = full path of the python script.
## $runCheckSnpsScript   = full path of check_snps.sh
## $splitFile            = split file used for parallel runs.
## $pileupFile           = pileup file created from the reference genome.
## $pileupExomeFile      = pileup file for exome.

echo "===== Setup other parameters ===== "
type samtools &>/dev/null  && isSamtools=1
if [ "$isSamtools" == 1 ]; then
  echo "===== Samtools found. samtools --version-only:`samtools --version-only`"
else
  echo "Error! Samtools not found"
  exit 1
fi

type bcftools &>/dev/null  && isBcftools=1
if [ "$isBcftools" == 1 ]; then
  echo "===== BCFtools found. bcftools --version-only:`bcftools --version-only`"
else
  echo "Error! BCFtools not found"
  exit 1
fi

type parallel &>/dev/null  && isParallel=1
if [ "$isParallel" == 1 ] && [ $parallelCount -gt 0 ]; then
  echo "===== Parallel mode: number of jobs=${parallelCount} ====="
  runCommand="parallel -j ${parallelCount}"
else
  echo "===== Non-Parallel mode: GNU parallel not found or invalid parallelCount ($parallelCount) ====="
  runCommand="xargs -I {}"
fi


if [ $numSplitFiles -le 0 ]; then
  echo "===== Number of split files (${numSplitFiles}) cannot be less than 0. Reset to 10 ====="
  numSplitFiles=10
fi



if [ -z ${dataDir} ] || [ ! -d ${dataDir} ]; then
  echo "ERROR! DataDir:${dataDir} does not exit. DatasetName=${datasetName}"
  exit 1
fi

if [ -z ${scriptDir} ] || [ ! -d ${scriptDir} ]; then
  echo "ERROR! scriptDir:${scriptDir} does not exit."
  exit 1
fi

if [ $isOriginalCaller -eq 1 ]; then
    echo "===== bcftools original mode ====="
    workingDir="${workingDir}/original/"
else
    echo "===== bcftools multiallec mode ====="
    workingDir="${workingDir}/multiallelic/"
fi

callsSingleDir="${workingDir}calls_single_vcfs/"
if [ ! -d ${callsSingleDir} ]; then
    mkdir -pv ${callsSingleDir}
fi

if [ $isTrio == 1 ];then
  pythonGetBaseCounts="${scriptDir}get_base_counts_trio.py"
  callsTrioDir="${workingDir}calls_trio_vcfs/"
  if [ ! -d ${callsTrioDir} ]; then
      mkdir -pv ${callsTrioDir}
  fi

else
  pythonGetBaseCounts="${scriptDir}get_base_counts_single.py"
fi

runCheckSnpsScript="${scriptDir}checkSnps.sh"
splitFile="${workindDir}split_chr${chromosome}_${datasetName}.split"
pileupFile="chr${chromosome}_${datasetName}.pileups"
pileupExomeFile="chr${chromosome}Ex_${datasetName}.pileups"



echo "===== Setup ====="
echo "DatasetName: ${datasetName}"
echo "WorkingDir: ${workingDir}"
echo "scriptDir: ${scriptDir}"
echo "dataDir: ${dataDir}"
echo "VariableSiteFile: ${variableSiteFile}"
echo "ReferenceGenome: ${refGenome}"
