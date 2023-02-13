rm(list=ls())

library(plyr)

# mr_results 
## adiposity_cancer ====
data <- fread("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_cancer"
adiposity_cancer <- data

data <- fread("analysis/004_colorectal_adiposity/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_cancer"
cancer_adiposity <- data

# adiposity_proteins ====
data <- fread("analysis/002_adiposity_proteins/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_proteins"
data <- data[,-c("sequence_id", "protein", "gene")]
adiposity_proteins <- data

data <- fread("analysis/005_proteins_adiposity/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_proteins"
data <- data[,-c("sequence_id", "protein", "gene")]
proteins_adiposity <- data

# proteins_cancer ====
data <- fread("analysis/003_proteins_colorectal/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
proteins_cancer <- data

data <- fread("analysis/006_colorectal_proteins/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
cancer_proteins <- data

# cis UVMR ====
data <- fread("analysis/011_protein_colorectal_cissnp/001_MR_results.txt", header = T)
data$analysis <- "cis-UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
cis <- data

# cis UVMR ====
data <- fread("analysis/010_UKB_proteins_colorectal/001_MR_results.txt", header = T)
data$analysis <- "UKB-UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
ukb <- data
head(cis)
head(data)

# master ====
data <- rbind.fill(adiposity_cancer, cancer_adiposity, adiposity_proteins, proteins_adiposity, proteins_cancer, cancer_proteins, cis, ukb)
write.table(data, "analysis/tables/UVMR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# mvmr_results 
## adiposity_cancer ====
data <- fread("analysis/008_mvmr/mvmr_results.txt", header = T)
data$analysis <- "MVMR"
data$analysis_ID <- "adiposity_protein_cancer"
data$mediator <- paste0(data$mediator, "_", data$protein, "_", data$gene)
mvmr <- data[,c("exposure", "mediator", "cancer", "group", "b", "se", "p", "OR", "lower_ci", "upper_ci", "Qstat", "Qpval", "fstat_adiposity", "fstat_protein")]
write.table(mvmr, "analysis/tables/MVMR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")









