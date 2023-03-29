rm(list=ls())
## set environment ====
library(dplyr)
library(readxl)

# data ====
data <- fread("analysis/010_UKB_proteins_colorectal/mr_results.txt")
data <- subset(data, outcome != "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt_with_header.txt")

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
data$group[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Female"
data$group[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "Male"

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
data$outcome[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "colon"
data$outcome[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "distal"
data$outcome[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "proximal"
data$outcome[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "rectal"
data$outcome[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "overall"
data$outcome[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt_with_header.txt"] <- "overall"

data$id.outcome <- paste0(data$outcome, "_", data$group)

# protein info ====
protein_info <- read_xlsx("../000_datasets/proteins/UKB_proteins_2022.xlsx", sheet = 7, skip = 2)
protein_info <- select(protein_info, `Assay Target`, `Target UniProt`)
protein_info <- unique(protein_info)
protein_info <- protein_info[!duplicated(protein_info$`Assay Target`),]
colnames(protein_info) <- c("exposure", "UNIPROT")
data <- left_join(data,protein_info, by = "exposure")

# methods ====
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame ====
data <- data[,c("exposure", "UNIPROT", "group", "outcome", "nsnp", "method", "b", "se", "pval", "OR", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/010_UKB_proteins_colorectal/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
