#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(maftools)
library(NMF)
library(MutationalPatterns)
library(dplyr)
library(readr)

RNA_maf_path <- args[1]
signatures <- args[2]
all_cont_plot <- args[3]
rel_cont_plot <- args[4]
file_path_96 <- args[5]

# The final output maf from STALLONE
RNA_maf <- read_tsv(RNA_maf_path)%>%
  filter(total_somatic_mutations > 15)

RNA_maf <- read.maf(RNA_maf)

# Read in the RNA maf as a TNM
RNA.tnm <- trinucleotideMatrix(maf = RNA_maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

RNA.mat <- RNA.tnm$nmf_matrix %>% t()

colnames(RNA.mat)<- substr(colnames(RNA.mat),1,12)

# Plot the signatures extracted from the TNM for each sample and save them
for (i in 1:ncol(RNA.mat)){
  
  sample <- colnames(RNA.mat)[i]
  sample_mat <- as.matrix(RNA.mat[,i])
  colnames(sample_mat) <- sample
  plot_96_profile(sample_mat)+
    ggsave(paste0(file_path_96, sample, ".pdf"))
  
}


# Read in the cosmic mutational signatures
cancer_signatures = read.table(args[2], sep = "\t", header = TRUE)

new_order = match(row.names(RNA.mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

# Fit the sample signatures to the cosmic signatures
fit_res_RNA <- fit_to_signatures(RNA.mat, cancer_signatures)

select <- which(rowSums(fit_res_RNA$contribution) > 0)

plot_contribution(fit_res_RNA$contribution[select,],
                  cancer_signatures[,select],
                  coord_flip = F,
                  mode = "absolute")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(all_cont_plot)

plot_contribution_heatmap(fit_res_RNA$contribution[select,],
                          cluster_samples = F,
                          method = "complete")+
  ggsave(rel_cont_plot)
