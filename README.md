# STALLONE
Somatic variant calling from a cohorts of unmatched RNA-Seq data. 

## Set up reference files. Run once.

### Set up the reference genome

### Make a STAR index
STAR --runMode genomeGenerate --genomeDir . \
--genomeFastaFiles GRCh37.fa --runThreadN 8

### Make a samtools index
samtools faidx GRCh37.fa

### Make a picard index
picard CreateSequenceDictionary R=GRCh37.fa O=human_g1k_v37.dict

### Create a chrom sizes file 
cut -f 1-2 GRCh37.fa.fai > GRCh37.fa.chromSizes

### Get the exons bed file ready
bedtools sort -i Homo_sapiens.GRCh37.87_exons_only.bed > exons1.bed
bedtools merge -i exons1.bed > exons_main_chrs.bed

## Description of scripts

VCF2MAF.sh converts the output vcfs into maf files
WES_vs_RNA_Seq.Rmd takes the raw GATK calls, compares them to WES from the TCGA and calcualtes mutational signatures/TMB
