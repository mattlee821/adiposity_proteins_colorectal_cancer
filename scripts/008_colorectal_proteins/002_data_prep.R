rm(list=ls())

# format results

# environment ====
## library ====
library(dplyr)
library(tidyr)

# data ====
data <- read.table("analysis/006_colorectal_proteins/mr_results.txt", header = T, sep = "\t")
data$lower_ci <- data$b - (1.96 * data$se)
data$upper_ci <- data$b + (1.96 * data$se)

# sort names out 
## groups 
table(data$exposure)
data$group[data$exposure == "overall_combined"] <- "Sex-combined"
data$group[data$exposure == "overall_combined_HRC"] <- "Sex-combined"
data$group[data$exposure == "colon_combined"] <- "Sex-combined"
data$group[data$exposure == "distal_combined"] <- "Sex-combined"
data$group[data$exposure == "proximal_combined"] <- "Sex-combined"
data$group[data$exposure == "rectal_combined"] <- "Sex-combined"

data$group[data$exposure == "overall_male"] <- "Male"
data$group[data$exposure == "colon_male"] <- "Male"
data$group[data$exposure == "distal_male"] <- "Male"
data$group[data$exposure == "proximal_male"] <- "Male"
data$group[data$exposure == "rectal_male"] <- "Male"

data$group[data$exposure == "overall_female"] <- "Female"
data$group[data$exposure == "colon_female"] <- "Female"
data$group[data$exposure == "distal_female"] <- "Female"
data$group[data$exposure == "proximal_female"] <- "Female"
data$group[data$exposure == "rectal_female"] <- "Female"

## exposures
table(data$exposure)
data$exposure[data$exposure == "overall_combined"] <- "overall"
data$exposure[data$exposure == "overall_combined_HRC"] <- "overall_HRC"
data$exposure[data$exposure == "colon_combined"] <- "colon"
data$exposure[data$exposure == "distal_combined"] <- "distal"
data$exposure[data$exposure == "proximal_combined"] <- "proximal"
data$exposure[data$exposure == "rectal_combined"] <- "rectal"

data$exposure[data$exposure == "overall_male"] <- "overall"
data$exposure[data$exposure == "colon_male"] <- "colon"
data$exposure[data$exposure == "distal_male"] <- "distal"
data$exposure[data$exposure == "proximal_male"] <- "proximal"
data$exposure[data$exposure == "rectal_male"] <- "rectal"

data$exposure[data$exposure == "overall_female"] <- "overall"
data$exposure[data$exposure == "colon_female"] <- "colon"
data$exposure[data$exposure == "distal_female"] <- "distal"
data$exposure[data$exposure == "proximal_female"] <- "proximal"
data$exposure[data$exposure == "rectal_female"] <- "rectal"

# format protein names ====
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
data <- data[, -which(names(data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10", "col11"))]

## methods
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame
data <- data[,c("exposure", "group", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b", "se", "pval", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/006_colorectal_proteins/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
