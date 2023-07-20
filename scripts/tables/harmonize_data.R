rm(list=ls())
# MR analysis of measures of adiposity and endometrial cancer 

# environment ====
## library ====
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(dplyr)
library(tidyverse)

# data ====
file_paths <- list.files(path = "analysis/", pattern = "harmonise", recursive = TRUE, full.names = TRUE)
file_paths <- file_paths[-grep("008_mvmr", file_paths)]
files <- lapply(file_paths, read.table, header = TRUE) # modify the read.table parameters as necessary

data <- do.call(bind_rows, files)

cols <- c("SNP","chr.exposure","pos.exposure","chr.outcome","pos.outcome",
          "effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","eaf.exposure","eaf.outcome",
          "beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome","samplesize.exposure","samplesize.outcome",
          "id.exposure","exposure","id.outcome","outcome",
          "palindromic","ambiguous","remove")
data <- data[, !names(data) %in% cols]

write.table(data, "analysis/tables/harmonize_data_all.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
