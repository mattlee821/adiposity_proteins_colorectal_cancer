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

### colours
#install.packages("wesanderson")
library(wesanderson)
d1 <- wes_palette("Royal1", type = "discrete")
d2 <- wes_palette("GrandBudapest2", type = "discrete")
d3 <- wes_palette("Cavalcanti1", type = "discrete")
d4 <- wes_palette("Rushmore1", type = "discrete")
discrete_wes_pal <- c(d1, d2, d3, d4)
rm(d1,d2,d3,d4)

## extract exposure instruments ====
data <- read_exposure_data("data/000_protein_data/exposure_data/cis_snps/cis_snps.txt",
                           clump = F,
                           sep = "\t",
                           snp_col = "rsID",
                           beta_col = "BETA",
                           se_col = "SE",
                           eaf_col = "EAF",
                           effect_allele_col = "EffectAllele",
                           other_allele_col = "OtherAllele",
                           pval_col = "Pval",
                           samplesize_col = "N",
                           phenotype = "phenotype",
                           min_pval = 5e-9)

data <- data[!duplicated(data[,"exposure"]),]

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
# remove extra cols 
data <- data[, -which(names(data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"))]

# exposure data ====
data$exposure <- as.factor(data$exposure)
data$f_stats <- (data$beta.exposure / data$se.exposure)^2 
exposure_data <- data
write.table(exposure_data, "analysis/011_protein_colorectal_cissnp/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# outcome data ====
## extract outcome data ====
filenames <- list.files(path = "../000_datasets/colorectal_cancer/GECCO_125K_GWAS_results/", pattern="annotated.txt", full.names=F)
list <- lapply(paste0("../000_datasets/colorectal_cancer/GECCO_125K_GWAS_results/",filenames), 
               read_outcome_data,
               snps = exposure_data$SNP,
               sep = " ",
               snp_col = "SNP",
               beta_col = "Effect",
               se_col = "StdErr",
               eaf_col = "Freq1",
               effect_allele_col = "Allele1",
               other_allele_col = "Allele2",
               pval_col = "P.value",
               min_pval = 1e-200,
               log_pval = FALSE,
               chr_col = "Chr",
               pos_col = "Position",
               phenotype_col = "phenotype"
)

# add outcome info
for (i in 1:length(filenames)){
  list[[i]]$outcome <- filenames[i]
}
# creat single dataframe
outcome_data <- bind_rows(list)
write.table(outcome_data, "analysis/011_protein_colorectal_cissnp/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
write.table(harmonise_data, "analysis/011_protein_colorectal_cissnp/harmonise_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
mr_results <- mr(harmonise_data)
write.table(mr_results, "analysis/011_protein_colorectal_cissnp/mr_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




