rm(list=ls())

# environment ====
library(dplyr)
library(data.table)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(tidyr)
library(ggforestplot)

# data ====
coloc <- fread("analysis/009_colocalisation/results/000_coloc_results.txt")
length(unique(coloc$exposure))
mr <- fread("analysis/011_protein_colorectal_cissnp/002_MR_results_shared_with_MVMR.txt")
length(unique(mr$exposure))

# format ====
# sequence_id column 
ncols <- max(stringr::str_count(coloc$exposure, "_")) + 1
colmn <- paste0("col", 1:ncols)
coloc <-
  tidyr::separate(
    data = coloc,
    col = exposure,
    sep = "_",
    into = colmn,
    remove = FALSE)
coloc$sequence_id <- paste0(coloc$col1, "_", coloc$col2)
coloc <- select(coloc, !colmn)

# shared ====
coloc_id <- unique(coloc$sequence_id)
mr_id <- unique(mr$sequence_id)
shared <- intersect(coloc_id, mr_id)
not_shared <- setdiff(coloc_id, mr_id)

a <- mr[mr$sequence_id %in% shared, ]
a <- coloc[coloc$sequence_id %in% shared, ]
