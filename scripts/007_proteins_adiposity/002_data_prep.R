rm(list=ls())
## set environment ====
library(dplyr)

# data ====
data <- read.table("analysis/005_proteins_adiposity/mr_results.txt", header = T, sep = "\t")

# CI
data$lower_ci <- data$b - (1.96 * data$se)
data$upper_ci <- data$b + (1.96 * data$se)

# format adiposity names ====
table(data$outcome)
data$group[data$outcome == "bmi_combined.txt"] <- "Sex-combined"
data$group[data$outcome == "whr_combined.txt"] <- "Sex-combined"

data$group[data$outcome == "bmi_male.txt"] <- "Male"
data$group[data$outcome == "whr_male.txt"] <- "Male"

data$group[data$outcome == "bmi_female.txt"] <- "Female"
data$group[data$outcome == "whr_female.txt"] <- "Female"

## outcome
data$outcome[data$outcome == "bmi_combined.txt"] <- "BMI"
data$outcome[data$outcome == "whr_combined.txt"] <- "WHR"

data$outcome[data$outcome == "bmi_male.txt"] <- "BMI"
data$outcome[data$outcome == "whr_male.txt"] <- "WHR"

data$outcome[data$outcome == "bmi_female.txt"] <- "BMI"
data$outcome[data$outcome == "whr_female.txt"] <- "WHR"

data$id.outcome <- paste0(data$outcome, "_", data$group)

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
data <- subset(data, select= -c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10))

# methods ====
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame ====
data <- data[,c("exposure", "group", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/005_proteins_adiposity/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
