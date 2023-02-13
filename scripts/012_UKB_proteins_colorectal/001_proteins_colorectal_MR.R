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
library(readxl)

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
exposure_data <- read_xlsx("../000_datasets/proteins/UKB_proteins_2022.xlsx", sheet = 7, skip = 2)
a <- as.data.frame(str_split(exposure_data$`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, pattern = ":", n = Inf, simplify = T))
colnames(a) <- c("CHR", "BP", "OA", "EA", "IMP", "v1")
a <- select(a, CHR, BP, OA, EA)
exposure_data <- cbind(exposure_data,a)
exposure_data <- select(exposure_data, rsID, CHR, `GENPOS (hg38)`, EA,OA, `Assay Target`, `A1FREQ (discovery)`, `BETA (discovery, wrt. A1)`, `SE (discovery)`, `log10(p) (discovery)`)
exposure_data$pval <- exp(exposure_data$`log10(p) (discovery)`)
colnames(exposure_data) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "exposure", "eaf.exposure", "beta.exposure", "se.exposure", "log10p.exposure", "p.exposure")
exposure_data$id.exposure <- exposure_data$exposure

exposure_data$exposure <- as.factor(exposure_data$exposure)
exposure_data$f_stats <- (exposure_data$b / exposure_data$se)^2 
exposure_data %>%
  group_by(exposure) %>%
  summarise(mean = mean(f_stats))
write.table(exposure_data, "analysis/010_UKB_proteins_colorectal/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

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
write.table(outcome_data, "analysis/010_UKB_proteins_colorectal/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data, action = 2)
write.table(harmonise_data, "analysis/010_UKB_proteins_colorectal/harmonise_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
mr_results <- mr(harmonise_data, method_list = methods)
write.table(mr_results, "analysis/010_UKB_proteins_colorectal/mr_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Sensitivity analysis ====
mr_singlesnp <- mr_singlesnp(harmonise_data)
write.table(mr_singlesnp, "analysis/010_UKB_proteins_colorectal/mr_singlesnp.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_hetrogeneity <- mr_heterogeneity(harmonise_data, method_list = methods_heterogeneity)
write.table(mr_hetrogeneity, "analysis/010_UKB_proteins_colorectal/mr_hetrogeneity.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_pleiotropy <- mr_pleiotropy_test(harmonise_data)
write.table(mr_pleiotropy, "analysis/010_UKB_proteins_colorectal/mr_pleiotropy.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

mr_leaveoneout <- mr_leaveoneout(harmonise_data)
write.table(mr_leaveoneout, "analysis/010_UKB_proteins_colorectal/mr_leaveoneout.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Plots ====
source("scripts/my_mr_scatter_plot.R")
plot_mr_scatter <- my_mr_scatter_plot(mr_results, harmonise_data)

plot_singlesnp_forest <- mr_forest_plot(mr_singlesnp)

plot_leaveoneout_forest <- mr_leaveoneout_plot(mr_leaveoneout)

plot_mr_funnel <- mr_funnel_plot(mr_singlesnp)

### save plots ====
pdf("analysis/010_UKB_proteins_colorectal/figures/plot_mr_scatter.pdf")
for (i in 1:length(plot_mr_scatter)) {
  print(plot_mr_scatter[[i]])
}
dev.off()

pdf("analysis/010_UKB_proteins_colorectal/figures/plot_singlesnp_forest.pdf")
for (i in 1:length(plot_singlesnp_forest)) {
  print(plot_singlesnp_forest[[i]])
}
dev.off()

pdf("analysis/010_UKB_proteins_colorectal/figures/plot_leaveoneout_forest.pdf")
for (i in 1:length(plot_leaveoneout_forest)) {
  print(plot_leaveoneout_forest[[i]])
}
dev.off()

pdf("analysis/010_UKB_proteins_colorectal/figures/plot_mr_funnel.pdf")
for (i in 1:length(plot_mr_funnel)) {
  print(plot_mr_funnel[[i]])
}
dev.off()
