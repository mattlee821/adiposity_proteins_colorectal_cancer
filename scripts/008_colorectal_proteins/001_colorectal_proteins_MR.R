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

### methods
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods_heterogeneity <- methods_heterogeneity[c(1,2,3,5)]
methods <- methods[c(3,6,10,13),1]

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
exposure_data <- read_exposure_data("data/exposure_data/000_clumped_crc.txt",
                                  clump = F,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "beta.exposure",
                                  se_col = "se.exposure",
                                  eaf_col = "eaf.exposure",
                                  effect_allele_col = "effect_allele.exposure",
                                  other_allele_col = "other_allele.exposure",
                                  pval_col = "pval.exposure",
                                  samplesize_col = "samplesize.exposure",
                                  phenotype = "id.exposure",
                                  min_pval = 5e-8)

exposure_data$exposure <- as.factor(exposure_data$exposure)
exposure_data$id.exposure <- exposure_data$exposure
exposure_data$f_stats <- (exposure_data$b / exposure_data$se)^2 
exposure_data %>%
group_by(exposure) %>%
summarise(mean = mean(f_stats))
write.table(exposure_data, "analysis/006_colorectal_proteins/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## extract outcome data ====
outcome_data <- read.table("~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/crc_proteins.txt", 
                           header = F, sep = "\t", fill = T, col.names = c("CHR", "POS", "SNPID", "SNP", 
                                                                           "EA", "OA", "beta.outcome",
                                                                           "pval.outcome", "log10pval", "se.outcome", "samplesize.outcome",
                                                                           "impMAF", "outcome",
                                                                           "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"))

outcome_data$outcome <- gsub(pattern = "/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", replacement = "", outcome_data$outcome)
outcome_data$outcome <- gsub(pattern = ".txt.gz.unzipped", replacement = "", outcome_data$outcome)
outcome_data$id.outcome <- outcome_data$outcome
outcome_data$mr_keep.outcome <- TRUE

# sequence_id column 
ncols <- max(stringr::str_count(outcome_data$outcome, "_")) + 1
colmn <- paste0("col", 1:ncols)
outcome_data <-
  tidyr::separate(
    data = outcome_data,
    col = outcome,
    sep = "_",
    into = colmn,
    remove = FALSE)
outcome_data$sequence_id <- paste0(outcome_data$col1, "_", outcome_data$col2)

# protein_name column 
outcome_data$id.outcome <- sub(".*?_", "", outcome_data$outcome) # remove first part of sequence_id from column
outcome_data$id.outcome <- sub(".*?_", "", outcome_data$id.outcome) # remove second part of sequence_id from column
outcome_data$protein <- sub(".*?_", "", outcome_data$id.outcome) # remove genename

# gene name column 
outcome_data$gene <- gsub("\\_.*","",outcome_data$id.outcome)

# remove extra cols 
outcome_data <- outcome_data[, -which(names(outcome_data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10","col11"))]
write.table(outcome_data, "analysis/006_colorectal_proteins/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
write.table(harmonise_data, "analysis/006_colorectal_proteins/harmonise_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
mr_results <- mr(harmonise_data, method_list = methods)
write.table(mr_results, "analysis/006_colorectal_proteins/mr_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Sensitivity analysis ====
mr_singlesnp <- mr_singlesnp(harmonise_data)
write.table(mr_singlesnp, "analysis/006_colorectal_proteins/mr_singlesnp.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_hetrogeneity <- mr_heterogeneity(harmonise_data, method_list = methods_heterogeneity)
write.table(mr_hetrogeneity, "analysis/006_colorectal_proteins/mr_hetrogeneity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_pleiotropy <- mr_pleiotropy_test(harmonise_data)
write.table(mr_pleiotropy, "analysis/006_colorectal_proteins/mr_pleiotropy.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_leaveoneout <- mr_leaveoneout(harmonise_data)
write.table(mr_leaveoneout, "analysis/006_colorectal_proteins/mr_leaveoneout.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Plots ====
source("scripts/my_mr_scatter_plot.R")
plot_mr_scatter <- my_mr_scatter_plot(mr_results, harmonise_data)

plot_singlesnp_forest <- mr_forest_plot(mr_singlesnp)

plot_leaveoneout_forest <- mr_leaveoneout_plot(mr_leaveoneout)

plot_mr_funnel <- mr_funnel_plot(mr_singlesnp)

### save plots ====
pdf("analysis/006_colorectal_proteins/figures/plot_mr_scatter.pdf")
for (i in 1:length(plot_mr_scatter)) {
  print(plot_mr_scatter[[i]])
}

pdf("analysis/006_colorectal_proteins/figures/plot_singlesnp_forest.pdf")
for (i in 1:length(plot_singlesnp_forest)) {
  print(plot_singlesnp_forest[[i]])
}
dev.off()

pdf("analysis/006_colorectal_proteins/figures/plot_leaveoneout_forest.pdf")
for (i in 1:length(plot_leaveoneout_forest)) {
  print(plot_leaveoneout_forest[[i]])
}
dev.off()

pdf("analysis/006_colorectal_proteins/figures/plot_mr_funnel.pdf")
for (i in 1:length(plot_mr_funnel)) {
  print(plot_mr_funnel[[i]])
}
dev.off()
