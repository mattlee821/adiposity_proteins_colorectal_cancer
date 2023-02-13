rm(list=ls())
## set environment ====
library(dplyr)

# data ====
data <- fread("analysis/002_adiposity_proteins/mr_results.txt")

# CI
data$lower_ci <- data$b - (1.96 * data$se)
data$upper_ci <- data$b + (1.96 * data$se)

# format adiposity names ====
table(data$exposure)
data$group[data$exposure == "BMI_combined"] <- "Sex-combined"
data$group[data$exposure == "WHR_combined"] <- "Sex-combined"
data$group[data$exposure == "BMI_female"] <- "Female"
data$group[data$exposure == "WHR_female"] <- "Female"
data$group[data$exposure == "BMI_male"] <- "Male"
data$group[data$exposure == "WHR_male"] <- "Male"

## exposures
data$exposure[data$exposure == "BMI_combined"] <- "BMI"
data$exposure[data$exposure == "WHR_combined"] <- "WHR"
data$exposure[data$exposure == "BMI_female"] <- "BMI"
data$exposure[data$exposure == "WHR_female"] <- "WHR"
data$exposure[data$exposure == "BMI_male"] <- "BMI"
data$exposure[data$exposure == "WHR_male"] <- "WHR"

data$id.exposure <- paste0(data$exposure, "_", data$group)

# format protein names ====
data$outcome <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/", "", data$outcome)

# sequence_id column 
ncols <- max(stringr::str_count(data$outcome, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = outcome,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

# protein_name column 
data$id.outcome <- sub(".*?_", "", data$outcome) # remove first part of sequence_id from column
data$id.outcome <- sub(".*?_", "", data$id.outcome) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.outcome) # remove genename

# gene name column 
data$gene <- gsub("\\_.*","",data$id.outcome)

# remove extra cols 
data <- subset(data, select= -c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10))

# methods ====
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame ====
data <- data[,c("exposure", "group", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/002_adiposity_proteins/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
