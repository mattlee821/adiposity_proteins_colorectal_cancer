rm(list=ls())

# environment ====
library(dplyr)
library(data.table)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(tidyr)
library(ggforestplot)

# mvmr proteins ====
associated_proteins <- fread("analysis/007_adiposity_proteins_colorectal/phenospd/004_adiposity_protein_cancer_pos_neg.txt")
associated_proteins <- unique(associated_proteins$sequence_id)

mvmr_results <- fread("analysis/008_mvmr/mvmr_results.txt")
a <- subset(mvmr_results, exposure == "BMI")
b <- subset(mvmr_results, exposure == "WHR")
mvmr_results <- rbind(a,b)
# sequence_id column 
ncols <- max(stringr::str_count(mvmr_results$ID, "_")) + 1
colmn <- paste0("col", 1:ncols)
mvmr_results <-
  tidyr::separate(
    data = mvmr_results,
    col = ID,
    sep = "_",
    into = colmn,
    remove = FALSE)
mvmr_results$sequence_id <- paste0(mvmr_results$col2, "_", mvmr_results$col3)
mvmr_results <- select(mvmr_results, !colmn)

mvmr_results <- mvmr_results[mvmr_results$sequence_id %in% associated_proteins,]
length(unique(mvmr_results$sequence_id))
mvmr_proteins <- unique(mvmr_results$mediator)

# data ====
mr_results_cis <- fread("analysis/011_protein_colorectal_cissnp/001_MR_results.txt")

# shared analyses
mr_results_cis$ID <- paste0(mr_results_cis$sequence_id, "_", mr_results_cis$outcome, "_", mr_results_cis$group)
mr_results_cis_ID <- unique(mr_results_cis$ID)
mvmr_results$ID <- paste0(mvmr_results$sequence_id, "_", mvmr_results$cancer, "_", mvmr_results$group)
mvmr_results_ID <- unique(mvmr_results$ID)
mr_mvmr_shared <- intersect(mr_results_cis_ID, mvmr_results_ID)

mr_results_cis <- mr_results_cis[mr_results_cis$ID %in% mr_mvmr_shared, ]
mvmr_results <- mvmr_results[mvmr_results$ID %in% mr_mvmr_shared, ]

length(unique(mr_results_cis$ID))
length(unique(mvmr_results$ID))

# join cis and non cis results ====
mr_results <- fread("analysis/003_proteins_colorectal/001_MR_results.txt")
mr_results$ID <- paste0(mr_results$sequence_id, "_", mr_results$outcome, "_", mr_results$group)
mr_results <- subset(mr_results, method == "IVW-MRE")
mr_results <- mr_results[mr_results$ID %in% mr_mvmr_shared, ]

length(unique(mr_results$exposure))

# save
write.table(mr_results, "analysis/011_protein_colorectal_cissnp/002_MR_results_shared_with_MVMR.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
