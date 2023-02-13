rm(list=ls())
## set environment ====
library(dplyr)

# data ====
data <- fread("analysis/011_protein_colorectal_cissnp/mr_results.txt")

# OR and CI
data$OR <- exp(data$b)
data$lower_ci <- exp(data$b - (1.96 * data$se))
data$upper_ci <- exp(data$b + (1.96 * data$se))

# protein names from other file ====
a <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")
a <- a[,c("sequence_id", "protein", "gene")]
a <- unique(a)

## groups 
table(data$outcome)
data$group[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Female"

data$group[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Male"

data$group[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "joint_Early_Onset_wald_MAC50_1.txt.annotated.txt"] <- "Sex-combined"
table(data$group)

# outcomes ====
table(data$outcome)
data$outcome[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$outcome[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$outcome[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$outcome[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$outcome[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$outcome[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$outcome[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$outcome[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$outcome[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "overall"
data$outcome[data$outcome == "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt_with_header.txt"] <- "overall-HRC-EUR"
data$outcome[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "colon"
data$outcome[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "distal"
data$outcome[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "proximal"
data$outcome[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "rectal"
data$outcome[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "overall"
data$outcome[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "overall"
data$outcome[data$outcome == "joint_Early_Onset_wald_MAC50_1.txt.annotated.txt"] <- "early_onset"
table(data$outcome)

data$id.outcome <- paste0(data$outcome, "_", data$group)

# format protein names ====
# sequence_id column 
ncols <- max(stringr::str_count(data$exposure, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = exposure,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

# protein_name column 
data$id.exposure <- sub(".*?_", "", data$exposure) # remove first part of sequence_id from column
data$id.exposure <- sub(".*?_", "", data$id.exposure) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.exposure) # remove genename

# gene name column 
data$gene <- gsub("\\_.*","",data$id.exposure)

# remove extra cols 
data <- select(data, !colmn)

# organise data frame ====
data <- data[,c("exposure", "group", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "OR", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/011_protein_colorectal_cissnp/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
