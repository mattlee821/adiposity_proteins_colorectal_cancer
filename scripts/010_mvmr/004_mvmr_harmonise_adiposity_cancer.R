rm(list=ls())
# MR analysis of measures of adiposity and metabolites 

# environment ====
## library ====
library(TwoSampleMR)
library(data.table)
library(dplyr)

# instruments ====
adiposity <- read.table("data/000_mvmr_data/001_adiposity_data.txt", header = T, sep = "\t")
cancer <- read.table("data/000_mvmr_data/001_cancer_data.txt", header = T, sep = "\t")

# harmonize adiposity:cancer
harmonise_data <- harmonise_data(adiposity, cancer, action = 2)
harmonise_data <- subset(harmonise_data, (as.character(id.exposure) == as.character(id.outcome))) # make sure its only the same ids matching because the harmonise will do ALL
write.table(harmonise_data, "analysis/008_mvmr/harmonise_data_adiposity_cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
