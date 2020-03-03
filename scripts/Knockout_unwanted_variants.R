#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)

# Parse the input args
TMB_out <- args[1]
maf_out <- args[2]
maf_files <- args[grepl("\\.maf",args)]
coverage_files <- args[grepl("_bases_covered.txt",args)]

# Read in all the mafs
bind <- list()
TMB <- list()
for (i in 1:length(maf_files)){
  
  file <- read_tsv(maf_files[[i]],col_names = T,comment = "#",col_types = cols(.default = "c"))%>%
    mutate(t_depth = as.numeric(t_depth))%>%
    mutate(t_alt_count = as.numeric(t_alt_count))%>%
    mutate(VAF = t_alt_count/t_depth)%>%
    mutate(Start_Position = as.numeric(Start_Position))%>%
    # QD < 2 is suggested by GATK to be a good cutoff
    mutate(QD = as.numeric(QD))%>%
    # FS > 30 is suggested by GATK to be a good cutoff
    mutate(FS = as.numeric(FS))
  
  bind[[i]] <- file
  
  cov <- read_tsv(coverage_files[[i]],col_names = F)
  colnames(cov) <- "Bases_sufficient_coverage"
  cov$Tumor_Sample_Barcode <- unique(file$Tumor_Sample_Barcode)
  
  coverage[[i]] <- cov
  
}
bound_coverage <- bind_rows(coverage)
# Mutation types to remove
remove <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "Silent")

# Bind together mafs
bound_mafs <- bind_rows(bind)%>%
  # Keep only protein coding variants
  filter(BIOTYPE == "protein_coding")%>%
  filter(!Variant_Classification %in% remove)%>%
  mutate(in_dbsnp = ifelse(dbSNP_RS != "novel", "In dbSNP", "Not in dbSNP"))%>%
  mutate(in_dbsnp = replace(in_dbsnp, is.na(in_dbsnp), "Not in dbSNP"))

# Read in the radar DB and set it up to be joined to the start position of the maf
radarDB <- read_tsv("Reference/RADAR_Human_AG_all_hg19_v2.txt")%>%
  select(chromosome, Start_Position = position)%>%
  mutate(Chromosome = gsub("chr", "", chromosome))%>%
  mutate(In_db = "Yes")

maf_filtered <- bound_mafs %>%
  # Join in the RADAR data
  left_join(radarDB)%>%
  mutate(in_RADAR = ifelse(!is.na(In_db), "In RADAR", "Not in RADAR"))%>%
  # Remove dbsnp reads
  filter(in_dbsnp == "Not in dbSNP")%>%
  # Remove anything in ExAc
  filter(is.na(ExAC_AF_Adj))%>%
  group_by(Chromosome, Start_Position)%>%
  mutate(repeated = n())%>%
  ungroup()%>%
  # Mutations repeated more than twice should be removed 
  filter(repeated <=2)%>%
  # Remove stuff in the RADAR db
  # The majority of RNA editing is in the UTRs so this won’t be 
  # as much of an issue if you stick with coding mutations only
  filter(in_RADAR == "Not in RADAR")%>%
  # Filter complete homozygous or perfect hets
  filter(VAF <1)%>%
  filter(VAF != 0.5)%>%
  filter(QD > 4)%>%
  filter(FS < 32)%>%
  # Remove low AF things
  filter(VAF > 0.1)%>%
  # Remove a variant allele count < 3
  filter(t_alt_count > 3)%>%
  # Frame shift insertions and splice region variants are most often wrong due to the nature of aligning RNA
  filter(! Variant_Classification %in% c("Frame_Shift_Ins", "Splice_Region"))%>%
  # Join in the coverage information
  left_join(bound_coverage)%>%
  group_by(Tumor_Sample_Barcode)%>%
  # Get the total number of mutations in each sample
  mutate(total_somatic_mutations= n())%>%
  # Calculate TMB
  mutate(TMB = (total_somatic_mutations/Bases_sufficient_coverage) * 1E6)%>%
  write_tsv(maf_out)










