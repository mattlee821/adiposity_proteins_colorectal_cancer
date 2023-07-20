rm(list=ls())

library(plyr)
library(dplyr)
library(data.table)
library(openxlsx)

# mr_results 
## adiposity_cancer ====
data <- fread("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T)
data <- subset(data, outcome != "overall-HRC-EUR")
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_cancer"
forward <- data

data <- fread("analysis/004_colorectal_adiposity/001_MR_results.txt", header = T)
data <- subset(data, exposure != "overall_HRC")
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_cancer"
reverse <- data

### join and save 
names(forward)
names(reverse)
forward <- select(forward, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(forward) <- c("adiposity", "cancer", "group", "nsnp_adiposity", "method", "b_forward", "se_forward", "pval_forward")
reverse <- select(reverse, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(reverse) <- c("cancer", "adiposity", "group", "nsnp_cancer", "method", "b_reverse", "se_reverse", "pval_reverse")
data <- left_join(forward,reverse)
head(data)
table_adiposity_cancer <- select(data, adiposity, cancer, group, nsnp_adiposity, nsnp_cancer, method, b_forward, se_forward, pval_forward, b_reverse, se_reverse, pval_reverse)

# adiposity_proteins ====
data <- fread("analysis/002_adiposity_proteins/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_proteins"
data <- data[,-c("sequence_id", "protein", "gene")]
forward <- data

data <- fread("analysis/005_proteins_adiposity/001_MR_results.txt", header = T)
data$analysis <- "UVMR"
data$analysis_ID <- "adiposity_proteins"
data <- data[,-c("sequence_id", "protein", "gene")]
reverse <- data

### join and save 
names(forward)
names(reverse)
forward <- select(forward, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(forward) <- c("adiposity", "protein", "group", "nsnp_adiposity", "method", "b_forward", "se_forward", "pval_forward")
reverse <- select(reverse, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(reverse) <- c("protein", "adiposity", "group", "nsnp_protein", "method", "b_reverse", "se_reverse", "pval_reverse")
data <- left_join(forward,reverse)
head(data)
table_adiposity_protein <- select(data, adiposity, protein, group, nsnp_adiposity, nsnp_protein, method, b_forward, se_forward, pval_forward, b_reverse, se_reverse, pval_reverse)

# protein_cancer ====
data <- fread("analysis/011_protein_colorectal_cissnp/001_MR_results.txt", header = T)
data <- subset(data, outcome != "overall-HRC-EUR")
data <- subset(data, outcome != "early_onset")
data$analysis <- "cis-UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
forward <- data

data <- fread("analysis/006_colorectal_proteins/001_MR_results.txt", header = T)
data <- subset(data, exposure != "overall-HRC-EUR")
data <- subset(data, exposure != "early_onset")
data$analysis <- "UVMR"
data$analysis_ID <- "proteins_cancer"
data <- data[,-c("sequence_id", "protein", "gene")]
reverse <- data

### join and save 
names(forward)
names(reverse)
forward <- select(forward, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(forward) <- c("protein", "cancer", "group", "nsnp_protein", "method", "b_forward", "se_forward", "pval_forward")
reverse <- select(reverse, exposure, outcome, group, nsnp, method, b, se, pval)
colnames(reverse) <- c("cancer", "protein", "group", "nsnp_cancer", "method", "b_reverse", "se_reverse", "pval_reverse")
data <- left_join(forward,reverse, by = c("protein", "cancer", "group"))
head(data)
table_protein_cancer <- select(data, protein, cancer, group, nsnp_protein, nsnp_cancer, method.y, b_forward, se_forward, pval_forward, b_reverse, se_reverse, pval_reverse)
colnames(table_protein_cancer)[colnames(table_protein_cancer) == "method.y"] <- "method_reverse"

### colocalization
data <- read.table("analysis/009_colocalisation/results/002_coloc_results.txt", header = T, sep = "\t")
table_colocalization <- select(data, -id)
colnames(table_colocalization) <- c("protein", "cancer", "group", "nsnp_colocalization", "h0", "h1", "h2", "h3", "h4")
write.table(table_colocalization, "analysis/tables/colocalization.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
table_protein_cancer_coloc <- left_join(table_protein_cancer, table_colocalization, by = c("protein", "cancer", "group"))

# mvmr_results 
## adiposity_cancer ====
data <- fread("analysis/008_mvmr/mvmr_results.txt", header = T)
data$analysis <- "MVMR"
data$analysis_ID <- "adiposity_protein_cancer"
data <- subset(data, protein == "GREM1")
data <- subset(data, mediator != "BMI")
mvmr <- select(data, adiposity, protein, cancer, group, b, se, p, Qstat, Qpval, fstat_adiposity, fstat_protein, analysis)
mvmr_adiposity_cancer_id <- paste0(mvmr$adiposity, "_", mvmr$cancer, "_", mvmr$group)
uvmr <- table_adiposity_cancer
uvmr$id <- paste0(uvmr$adiposity, "_", uvmr$cancer, "_", uvmr$group)
uvmr <- uvmr[uvmr$id %in% mvmr_adiposity_cancer_id, ]
uvmr <- subset(uvmr, method == "IVW-MRE")
uvmr <- select(uvmr, adiposity, cancer, group, method, b_forward, se_forward, pval_forward)
colnames(uvmr) <- c("adiposity", "cancer", "group", "method_UVMR", "b_UVMR", "se_UVMR", "pval_UVMR")
table_mvmr_uvmr <- left_join(mvmr, uvmr, by = c("adiposity", "cancer", "group"))
table_mvmr_uvmr <- table_mvmr_uvmr[order(table_mvmr_uvmr$adiposity, 
                                         table_mvmr_uvmr$protein, 
                                         table_mvmr_uvmr$cancer, 
                                         table_mvmr_uvmr$group), ]


## combine ====
list_of_datasets <- list("UVMR; adiposity_cancer" = table_adiposity_cancer, 
                         "UVMR; adiposity_protein" = table_adiposity_protein,
                         "UVMR; protein_cancer" = table_protein_cancer_coloc,
                         "colocalization" = table_colocalization,
                         "MVMR" = table_mvmr_uvmr)
openxlsx::write.xlsx(list_of_datasets, file = "analysis/tables/results.xlsx")




