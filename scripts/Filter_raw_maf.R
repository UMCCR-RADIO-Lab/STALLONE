#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(readr)

maf_file <- args[1]
output_maf <- args[2]

# Read in the maf file
unfiltered_maf <- read_tsv(maf_file,col_names = T,comment = "#",col_types = cols(.default = "c"))

# Filter anything that is obviously a common varaint
filtered_maf <- unfiltered_maf %>%
    filter(FILTER != "common_variant")%>%
    write_tsv(output_maf)
