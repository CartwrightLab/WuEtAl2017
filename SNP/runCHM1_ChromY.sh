#!/usr/bin/bash
# This is the setup script of WuEtAl's analyses
# variables and directories need to change accordingly
#
# This will take one bam file,and generate a list of sites for which the kid is heterozygous
#
#
# Location of data file:
# CHM1 - run:SRR1283824, samples:SRS605680, experiment:SRX540665, study:SRP017546
# https://www.ncbi.nlm.nih.gov/sra/SRX540665[accn]
# ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283824/SRR1283824.sra
#
# Usage:
# runCHM1_ChromY 21


datasetName="CHM1"
chromosome=$1  #21
## sub chrome 10
if [ $chromosome -eq 10 ]; then
  chromosomePosition=":85534747-135534747"
fi

workingDir="./CHM1/${chm1Name}_C${chromosome}/"
scriptDir="./"
variableSiteFile="/storage/CEU/ALL.chr${chromosome}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
refGenome="/storage/reference_genomes/human/human_g1k_v37.fasta"

parallelCount=20
numSplitFiles=500
isOriginalCaller=1 #[1|0] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
isTrio=0
exome=0
exome_file=0

dataDir="/storage/chm1/bam/"
dataFiles="SRR1283824.bam"

# RUN
. ${scriptDir}pipelineCore.sh
