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
data <- read.table("analysis/007_adiposity_proteins_colorectal/harmonize_data_all.txt", header = T, sep = "\t")
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods_heterogeneity <- methods_heterogeneity[c(1,2,3,5)]
methods <- methods[c(3,6,10,13),1]

## Sensitivity analysis ====
mr_singlesnp <- mr_singlesnp(data)
write.table(mr_singlesnp, "analysis/007_adiposity_proteins_colorectal/mr_singlesnp.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
rm(mr_singlesnp)

mr_pleiotropy <- mr_pleiotropy_test(data)
write.table(mr_pleiotropy, "analysis/007_adiposity_proteins_colorectal/mr_pleiotropy.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
rm(mr_pleiotropy)

mr_leaveoneout <- mr_leaveoneout(data)
write.table(mr_leaveoneout, "mr_leaveoneout.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
rm(mr_leaveoneout)


## remove this pair and run
#'Bm0HTx' on 'yUtqQo'

mr_hetrogeneity <- mr_heterogeneity(data, method_list = methods_heterogeneity)
write.table(mr_hetrogeneity, "analysis/007_adiposity_proteins_colorectal/mr_hetrogeneity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
rm(mr_hetrogeneity)
