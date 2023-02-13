# clump instruments
rm(list=ls())

# environment ====
## library ====
# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
library(dplyr)

# read data
filenames <- list.files(path = "data/exposure_data/", pattern="*snps", full.names=F)
list <- lapply(paste0("data/exposure_data/",filenames), read.table, header = T)
list1 <- list[c(1:8)]
list2 <- list[c(9:16)]

# format
data <- bind_rows(list1, .id = "column_label")
data <- data[,-1]
data <- data[,c("SNP", "MarkerName", "chrom_start", "chr_name", "Allele1", "Allele2", "Freq1", "FreqSE",
 "MinFreq", "MaxFreq", "Effect", "StdErr", "P.value",
 "phenotype")]

data2 <- bind_rows(list2, .id = "column_label")
data2 <- data2[,-1]
data2 <- data2[,c("SNP", "MarkerName", "chrom_start", "chr_name", "Allele1", "Allele2", "Freq1", "FreqSE",
 "MinFreq", "MaxFreq", "Effect", "StdErr", "P.value",
 "phenotype")]

data <- bind_rows(data,data2)

## exposures
colnames(data)[5] <- "effect_allele.exposure"
colnames(data)[6] <- "other_allele.exposure"
colnames(data)[7] <- "eaf.exposure"
colnames(data)[11] <- "beta.exposure"
colnames(data)[12] <- "se.exposure"
colnames(data)[13] <- "pval.exposure"
colnames(data)[14] <- "exposure"
data$id.exposure <- data$exposure

data$exposure[data$exposure == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- "overall"
data$exposure[data$exposure == "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt"] <- "overall_HRC"
data$exposure[data$exposure == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$exposure[data$exposure == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$exposure[data$exposure == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$exposure[data$exposure == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$exposure[data$exposure == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$exposure[data$exposure == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$exposure[data$exposure == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$exposure[data$exposure == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$exposure[data$exposure == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- "colon"
data$exposure[data$exposure == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- "distal"
data$exposure[data$exposure == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- "overall"
data$exposure[data$exposure == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- "overall"
data$exposure[data$exposure == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- "proximal"
data$exposure[data$exposure == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- "rectal"

data$id.exposure[data$id.exposure == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- "overall_combined"
data$id.exposure[data$id.exposure == "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt"] <- "overall_combined_HRC"
data$id.exposure[data$id.exposure == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "colon_female"
data$id.exposure[data$id.exposure == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "colon_male"
data$id.exposure[data$id.exposure == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "distal_female"
data$id.exposure[data$id.exposure == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "distal_male"
data$id.exposure[data$id.exposure == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "proximal_female"
data$id.exposure[data$id.exposure == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "proximal_male"
data$id.exposure[data$id.exposure == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "rectal_female"
data$id.exposure[data$id.exposure == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "rectal_male"
data$id.exposure[data$id.exposure == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- "colon_combined"
data$id.exposure[data$id.exposure == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- "distal_combined"
data$id.exposure[data$id.exposure == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- "overall_female"
data$id.exposure[data$id.exposure == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- "overall_male"
data$id.exposure[data$id.exposure == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- "proximal_combined"
data$id.exposure[data$id.exposure == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- "rectal_combined"

write.table(data, "data/exposure_data/000_unclumped_crc.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- clump_data(data, clump_kb = 10000, clump_r2 = 0.001)

write.table(data, "data/exposure_data/000_clumped_crc.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
