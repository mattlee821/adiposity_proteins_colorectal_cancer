rm(list=ls())

# format results

# environment ====
## library ====
library(dplyr)
library(tidyr)

# data ====
filenames <- list.files(path = "analysis/001_adiposity_colorectal/", pattern="mr_results", full.names=F)
list <- lapply(paste0("analysis/001_adiposity_colorectal/",filenames), read.table, header = T, sep = "\t")
data <- bind_rows(list)

# OR and CI
data$OR <- exp(data$b)
data$lower_ci <- exp(data$b - (1.96 * data$se))
data$upper_ci <- exp(data$b + (1.96 * data$se))

# sort names out 
## groups 
table(data$exposure)
data$group[data$exposure == "BMI_combined"] <- "Sex-combined"
data$group[data$exposure == "WHR_combined"] <- "Sex-combined"

data$group[data$exposure == "BMI_male"] <- "Male"
data$group[data$exposure == "WHR_male"] <- "Male"

data$group[data$exposure == "BMI_female"] <- "Female"
data$group[data$exposure == "WHR_female"] <- "Female"

## exposures
data$exposure[data$exposure == "BMI_combined"] <- "BMI"
data$exposure[data$exposure == "WHR_combined"] <- "WHR"

data$exposure[data$exposure == "BMI_male"] <- "BMI"
data$exposure[data$exposure == "WHR_male"] <- "WHR"

data$exposure[data$exposure == "BMI_female"] <- "BMI"
data$exposure[data$exposure == "WHR_female"] <- "WHR"

## outcomes
table(data$outcome)
data$outcome[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$outcome[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$outcome[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$outcome[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$outcome[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$outcome[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$outcome[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$outcome[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$outcome[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- "overall"
data$outcome[data$outcome == "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt"] <- "overall-HRC-EUR"
data$outcome[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- "colon"
data$outcome[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- "distal"
data$outcome[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- "proximal"
data$outcome[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- "rectal"
data$outcome[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- "overall"
data$outcome[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- "overall"

## methods
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame
data <- data[,c("exposure", "outcome", "group", "nsnp", "method", "b", "se", "pval", "OR", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/001_adiposity_colorectal/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
