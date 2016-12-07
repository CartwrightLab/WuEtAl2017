#!/usr/bin/bash
# This is the setup script of WuEtAl's analyses
# variables and directories need to change accordingly
#
# This will take three bam files (child, dad, mom) and
# generate a list of sites for which the child is heterozygous
#
#
# Location of all files:
# Human Reference Genome:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
# Variable Site:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr10.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
# CEU2013:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12892/high_coverage_alignment/NA12892.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12891/high_coverage_alignment/NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
# CEU2012:
# ftp://ftp/technical/working/20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam
# ftp://ftp/technical/working/20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.20120117.bam
# ftp://ftp/technical/working/20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.20120117.bam
#
# CEU2011:
# ftp://ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam
# ftp://ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam
# ftp://ftp/technical/working/20110915_CEUtrio_b37_decoy_alignment/CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam
#
# CEU2010:
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100311.bam
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12892/alignment/NA12892.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12891/alignment/NA12891.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12878/alignment/NA12878.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100311.bam
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12892/alignment/NA12892.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam
# ftp://ftp/technical/pilot2_high_cov_GRCh37_bams/data/NA12891/alignment/NA12891.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam


#
# Usage:
# runCEU1X_ChromY.sh 13 21

datasetName="CEU${1}" #"CEU13"
chromosome=$2 #21
## sub chrome 10
if [ $chromosome -eq 10 ]; then
  chromosomePosition=":85534747-135534747"
fi

workingDir="./CEU/${datasetName}_C${chromosome}/"
scriptDir="./"
variableSiteFile="/storage/CEU/ALL.chr${chromosome}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
refGenome="/storage/reference_genomes/human/human_g1k_v37.fasta"

parallelCount=20
numSplitFiles=500
isOriginalCaller=1 #[1|0] # bcftool call:classic mode -c, --consensus-caller  or new mode -m, --multiallelic-caller
isTrio=1
exome=0
exome_file=0


if [ ${datasetName} == "CEU10" ]; then
  dataDir="/storage/CEUTrio/20100311/"
  dataFiles=('NA12878.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100311.bam' 'NA12892.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam' 'NA12891.chrom21.ILLUMINA.bwa.CEU.high_coverage.20100517.bam')
  if [ $chromosome -eq 10 ]; then
    dataFiles=('NA12878.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100311.bam' 'NA12891.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam' 'NA12892.chrom10.ILLUMINA.bwa.CEU.high_coverage.20100517.bam')
  fi

elif [ ${datasetName} == "CEU11" ]; then
  dataDir="/storage/CEUTrio/20110915/"
  dataFiles=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.bam')
elif [ ${datasetName} == "CEU12" ]; then
  dataDir="/storage/CEUTrio/20120117/"
  dataFiles=('CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12892.clean.dedup.recal.20120117.bam' 'CEUTrio.HiSeq.WGS.b37_decoy.NA12891.clean.dedup.recal.20120117.bam')
elif [ ${datasetName} == "CEU13" ]; then
  dataDir="/storage/CEUTrio/20130906/"
  dataFiles=('NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12892.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam' 'NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam')
fi


# RUN
. ${scriptDir}pipelineCore.sh
