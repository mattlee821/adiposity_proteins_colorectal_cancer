# make data frames for coloc
rm(list=ls())

# environment ====
## library ====
library(dplyr)
library(data.table)

# adiposity
data <- read.table("data/exposure_data/000_clumped_adiposity.txt", header = T, sep = "\t")
snps <- as.data.frame(unique(data$SNP))
dim(snps)

write.table(snps, "data/exposure_data/snp_list_adiposity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# crc
data <- read.table("data/exposure_data/000_clumped_crc.txt", header = T, sep = "\t")
snps <- as.data.frame(unique(data$SNP))
dim(snps)

write.table(snps, "data/exposure_data/snp_list_crc.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# proteins
data <- fread("data/exposure_data/000_clumped_proteins.txt", header = T, sep = "\t")
snps <- as.data.frame(unique(data$SNP))
dim(snps)

write.table(snps, "data/exposure_data/snp_list_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

