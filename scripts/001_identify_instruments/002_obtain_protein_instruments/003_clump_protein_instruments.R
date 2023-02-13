# clump instruments
rm(list=ls())

# environment ====
## library ====
# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(data.table)
setDTthreads(percent = 90)

# data
# data <- fread("data/exposure_data/protein_instruments.txt ", header = T, fill = T)
data <- read.table("data/exposure_data/protein_instruments.txt", header = F, fill = T,
  col.names = c("CHR", "POS", "SNPID", "SNP", 
  "EA", "OA", "beta.exposure", "pval.exposure", "minus_log10_pval", "se.exposure", "samplesize.exposure", "ImpMAF", 
  "exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure"))

data$exposure <- gsub(pattern = "/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", replacement = "", data$exposure)
data$exposure <- gsub(pattern = ".txt.gz.unzipped", replacement = "", data$exposure)
data$id.exposure <- data$exposure
write.table(data, "data/exposure_data/000_unclumped_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
length(unique(data$id.exposure))

# clump SNPs ====
data1 <- clump_data(data, clump_kb = 10000, clump_r2 = 0.001)
length(unique(data1$id.exposure))

write.table(data1, "data/exposure_data/000_clumped_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




