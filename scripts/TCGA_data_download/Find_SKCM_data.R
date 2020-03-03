library(tibble)
library(tidyr)
library(readr)
library(dplyr)
library(TCGAbiolinks)

# This script uses TCGA biolinks to get a random sample of 26 tumours that also have matched blood
# and WES information

# Nicen this script up if you go any further.

# Find the LUSC primary RNA-Seq
query_rseq <- GDCquery(project = c("TCGA-SKCM"),
                       data.category = "Sequencing Reads",  
                       sample.type = "Primary solid Tumor",experimental.strategy = "RNA-Seq")

# Get results and keep only primary samples
results_rseq <- getResults(query_rseq)%>%
  filter(!duplicated(substr(cases,1,12)))


rseq_manifest <- getManifest(query_rseq)

# Save a manifest of all 30 samples for me to download and test RSVC on
set.seed(100)
samples <- sample(2:nrow(rseq_manifest), 30)
to_save_30 <- rseq_manifest [samples,]

write_tsv(to_save_30, "/data/cephfs/punim0648/Pattison_projects/STALLONE/TCGA_manifests/SKCM_RNA_manifest.tsv",col_names = T)

# Find the LUSC primary RNA-Seq
query_rseq_normal <- GDCquery(project = c("TCGA-SKCM"),
                              data.category = "Sequencing Reads",  
                              sample.type = "Solid Tissue Normal",experimental.strategy = "RNA-Seq")
# Get results and keep only primary samples
results_rseq_normal <- getResults(query_rseq_normal)%>%
  filter(!duplicated(substr(cases,1,12)))

# Get the manifest of the full panel of normals for PON creation
rseq_normal_manifest_all <- getManifest(query_rseq_normal,save = FALSE)%>%
  filter(id %in% results_rseq_normal$id)

# Arrange so that the samples with matched adjacent lung are incorporated. 

results_rseq <- filter(results_rseq, substr(cases,1,12) %in% substr(results_rseq_normal$cases,1,12))

# Keep only the matched normal samples as well
results_rseq_normal <- results_rseq_normal %>%
  filter(substr(cases,1,12) %in% substr(results_rseq$cases,1,12))

# Get clinical data from TCGA biolinks
clinical <- GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")

# Filter clinical data for the patients 
clinical_filtered <- filter(clinical, submitter_id %in% substr(results_rseq$cases,1,12))
