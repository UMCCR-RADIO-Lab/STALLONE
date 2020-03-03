library(tibble)
library(tidyr)
library(readr)
library(dplyr)
library(TCGAbiolinks)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
# This script uses TCGA biolinks to get a random sample of 26 tumours that also have matched blood
# and WES information

# Nicen this script up if you go any further.

# Find the LUSC primary RNA-Seq
query_rseq <- GDCquery(project = c("TCGA-LUSC"),
                       data.category = "Sequencing Reads",  
                       sample.type = "Primary solid Tumor",experimental.strategy = "RNA-Seq")
# Get results and keep only primary samples
results_rseq <- getResults(query_rseq)%>%
  filter(!duplicated(substr(cases,1,12)))

# Save a manifest of 150 samples for me to download and test RSVC on
rseq_manifest <- getManifest(query_rseq)

set.seed(100)
samples <- sample(2:nrow(rseq_manifest), 150)

to_save_150 <- rseq_manifest [samples,]

write_tsv(to_save_150 , "/Users/adpattison/igv/mounts/Spartan_2/download-manifests/Tumour_RNA_Seq_manifest_150.tsv",col_names = T)

# Find the LUSC primary RNA-Seq
query_rseq_normal <- GDCquery(project = c("TCGA-LUSC"),
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

# Get the blood derived normals
#Blood normals

query_wxs <- GDCquery(project = c("TCGA-LUSC"),
                      data.category = "Sequencing Reads",
                      experimental.strategy = "WXS")

results_all_wxs <- getResults(query_wxs)

# This shows that there are only primary tumours, blood samples and matched tissue in this dataset

# barcode <- results_all_wxs %>%
#   select(cases)%>%
#   mutate(original = cases)
# 
# barcode_annotation_table <- separate(barcode, cases, into = c("Project", "Tissue_source_site", "Patient", "Sample_vial", "Portion_analyte", "Plate", "Center"),sep = "-")
# 
# table(barcode_annotation_table$Sample_vial)
# 
# table(barcode_annotation_table$Tissue_source_site)

# Get the WXS blood and tumour samples
bloods <- results_all_wxs%>%
  filter(substr(cases,14,15) == "10")%>%
  filter(!duplicated(substr(cases,1,12)))

# Get the manifest of the full panel of normals for PON creation
# 100 randomly selected normals
set.seed(100)
WXS_normal_manifest_PON <- getManifest(query_wxs,save = FALSE)%>%
  filter(id %in% bloods[sample(1:nrow(bloods),100),]$id)

tumours <-  results_all_wxs%>%
  filter(substr(cases,14,15) == "01")%>%
  filter(!duplicated(substr(cases,1,12)))

# I only want samples that have matched blood WXS as well
WXS_test <- tumours[substr(tumours$cases,1,12) %in% substr(bloods$cases,1,12),]
WXS_test <- WXS_test[substr(WXS_test$cases,1,12) %in% substr(results_rseq$cases,1,12),]

# Make sure we have matched samples down to the sample portion
results_wxs <- filter(WXS_test,substr(cases,1,15) %in% substr(results_rseq$cases,1,15))

bloods_keep <- bloods[substr(bloods$cases,1,12) %in% substr(results_wxs$cases,1,12),]

results_rseq <- filter(results_rseq,substr(cases,1,15) %in% substr(results_wxs$cases,1,15))

# Make sure have the same 26 samples in the normal manifest. 
results_rseq_normal <- filter(results_rseq_normal,substr(cases,1,12) %in% substr(results_wxs$cases,1,12))

# Check the samples are all the same
substr(bloods_keep$cases,1,12) %in% substr(results_rseq$cases,1,12)

wes_blood_manifest <- getManifest(query_wxs,save = FALSE) %>%
  filter(id %in% bloods_keep$id)

wes_manifest <- getManifest(query_wxs,save = FALSE) %>%
  filter(id %in% results_wxs$id)

rseq_manifest <- getManifest(query_rseq,save = FALSE)%>%
  filter(id %in% results_rseq$id)

rseq_normal_manifest <- getManifest(query_rseq_normal,save = FALSE)%>%
  filter(id %in% results_rseq_normal$id)

# Save the resulting matched manifest files 
write_tsv(wes_blood_manifest,"/Users/adpattison/igv/mounts/Spartan_2/manifests/WES_blood_manifest.tsv")
write_tsv(wes_manifest,"/Users/adpattison/igv/mounts/Spartan_2/manifests/WES_manifest.tsv")
write_tsv(rseq_manifest,"/Users/adpattison/igv/mounts/Spartan_2/manifests/RNA_Seq_manifest.tsv")
write_tsv(rseq_normal_manifest,"/Users/adpattison/igv/mounts/Spartan_2/manifests/RNA_Seq_normal_manifest.tsv")
write_tsv(rseq_normal_manifest_all,"/Users/adpattison/igv/mounts/Spartan_2/manifests/RNA_Seq_normal_manifest_all_samples.tsv")
write_tsv(WXS_normal_manifest_PON,"/Users/adpattison/igv/mounts/Spartan_2/manifests/WXS_normal_manifest_PON_samples.tsv")


# Get clinical data from TCGA biolinks
clinical <- GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")

# Filter clinical data for the patients 

clinical_filtered <- filter(clinical, submitter_id %in% substr(bloods_keep$cases,1,12))
