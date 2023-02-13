rm(list=ls())

# format results

# environment ====
## library ====
library(dplyr)
library(tidyr)

# data ====
filenames <- list.files(path = "analysis/004_colorectal_adiposity/", pattern="mr_results", full.names=F)
list <- lapply(paste0("analysis/004_colorectal_adiposity/",filenames), read.table, header = T, sep = "\t")
data <- bind_rows(list)

data$lower_ci <- data$b - (1.96 * data$se)
data$upper_ci <- data$b + (1.96 * data$se)

# sort names out 
## groups 
table(data$outcome)
data$group[data$outcome == "bmi_combined.txt"] <- "Sex-combined"
data$group[data$outcome == "whr_combined.txt"] <- "Sex-combined"

data$group[data$outcome == "bmi_male.txt"] <- "Male"
data$group[data$outcome == "whr_male.txt"] <- "Male"

data$group[data$outcome == "bmi_female.txt"] <- "Female"
data$group[data$outcome == "whr_female.txt"] <- "Female"

## outcomes
data$outcome[data$outcome == "bmi_combined.txt"] <- "BMI"
data$outcome[data$outcome == "whr_combined.txt"] <- "WHR"

data$outcome[data$outcome == "bmi_male.txt"] <- "BMI"
data$outcome[data$outcome == "whr_male.txt"] <- "WHR"

data$outcome[data$outcome == "bmi_female.txt"] <- "BMI"
data$outcome[data$outcome == "whr_female.txt"] <- "WHR"

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

## methods
data$method[data$method == "Inverse variance weighted (multiplicative random effects)"] <- "IVW-MRE"

# organise data frame
data <- data[,c("exposure", "outcome", "group", "nsnp", "method", "b", "se", "pval", "lower_ci", "upper_ci")]

# save
write.table(data, "analysis/004_colorectal_adiposity/001_MR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
