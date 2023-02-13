# clump instruments
rm(list=ls())

# environment ====
## library ====
# remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

# read data
filenames <- list.files(path = "data/exposure_data/", pattern="txt_pulit_snps.txt", full.names=F)
list <- lapply(paste0("data/exposure_data/",filenames), read.table, header = T)

# format and save ====
for (i in 1:length(filenames)){
  list[[i]]$adiposity <- filenames[i]
}
data <- bind_rows(list)
colnames(data) <- c("CHR","POS","SNP",
  "effect_allele.exposure","other_allele.exposure","eaf.exposure", 
  "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
  "INFO", "exposure", "id.exposure")

## exposures
data$exposure[data$exposure == "bmi_combined.txt"] <- "BMI"
data$exposure[data$exposure == "bmi_male.txt"] <- "BMI"
data$exposure[data$exposure == "bmi_female.txt"] <- "BMI"
data$exposure[data$exposure == "whr_combined.txt"] <- "WHR"
data$exposure[data$exposure == "whr_male.txt"] <- "WHR"
data$exposure[data$exposure == "whr_female.txt"] <- "WHR"

data$id.exposure[data$id.exposure == "bmi_combined.txt_pulit_snps.txt"] <- "BMI_combined"
data$id.exposure[data$id.exposure == "bmi_male.txt_pulit_snps.txt"] <- "BMI_male"
data$id.exposure[data$id.exposure == "bmi_female.txt_pulit_snps.txt"] <- "BMI_female"
data$id.exposure[data$id.exposure == "whr_combined.txt_pulit_snps.txt"] <- "WHR_combined"
data$id.exposure[data$id.exposure == "whr_male.txt_pulit_snps.txt"] <- "WHR_male"
data$id.exposure[data$id.exposure == "whr_female.txt_pulit_snps.txt"] <- "WHR_female"

write.table(data, "data/exposure_data/000_unclumped_adiposity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- clump_data(data, clump_kb = 10000, clump_r2 = 0.001)

write.table(data, "data/exposure_data/000_clumped_adiposity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
