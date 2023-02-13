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
exposure_data <- read_exposure_data("data/exposure_data/000_clumped_adiposity.txt",
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
                                  min_pval = 5e-9)

exposure_data$id.exposure <- exposure_data$exposure
exposure_data$exposure <- as.factor(exposure_data$exposure)
exposure_data$f_stats <- (exposure_data$b / exposure_data$se)^2 
exposure_data %>%
group_by(exposure) %>%
summarise(mean = mean(f_stats))

## extract outcome data ====
filenames <- list.files(path = "/data/GWAS_data/files/huyghe_2018_PMID30510241/processed/", pattern="annotated.txt", full.names=F)
list <- lapply(paste0("/data/GWAS_data/files/huyghe_2018_PMID30510241/processed/",filenames), 
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
outcome_data_combined <- bind_rows(list)

## harmonize data ====
harmonise_data <- harmonise_data(exposure_data, outcome_data_combined, action = 2)

# seperate out male, female, and combined into individual dataframes ====
male_exposure <- c("BMI_male","WHR_male")
female_exposure <- c("BMI_female","WHR_female")
combined_exposure <- c("BMI_combined","WHR_combined")

male_outcome <- c("joint_colon_Male_wald_MAC50_1.TBL.annotated.txt",
  "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt",
  "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt",
  "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt",
  "Stratified_male_Meta_125K_Results.tsv.annotated.txt")
female_outcome <- c("joint_colon_Female_wald_MAC50_1.TBL.annotated.txt",
  "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt",
  "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt",
  "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt",
  "Stratified_female_Meta_125K_Results.tsv.annotated.txt")
combined_outcome <- c("MarginalMeta_125K_Results.tsv.annotated.txt",
  "MarginalMeta_HRC_EUR_only_Results.tsv.annotated.txt",
  "Stratified_colon_Meta_125K_Results.tsv.annotated.txt",
  "Stratified_distal_Meta_125K_Results.tsv.annotated.txt",
  "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt",
  "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt",
  "joint_Early_Onset_wald_MAC50_1.txt.annotated.txt")

harmonise_data_male <- harmonise_data[harmonise_data$exposure %in% male_exposure, ]
harmonise_data_male <- harmonise_data_male[harmonise_data_male$outcome %in% male_outcome, ]

harmonise_data_female <- harmonise_data[harmonise_data$exposure %in% female_exposure, ]
harmonise_data_female <- harmonise_data_female[harmonise_data_female$outcome %in% female_outcome, ]

harmonise_data_combined <- harmonise_data[harmonise_data$exposure %in% combined_exposure, ]
harmonise_data_combined <- harmonise_data_combined[harmonise_data_combined$outcome %in% combined_outcome, ]

## MR ====
mr_results_male <- mr(harmonise_data_male, method_list = methods)
mr_results_female <- mr(harmonise_data_female, method_list = methods)
mr_results_combined <- mr(harmonise_data_combined, method_list = methods)

## Sensitivity analysis ====
mr_singlesnp_male <- mr_singlesnp(harmonise_data_male)
mr_singlesnp_female <- mr_singlesnp(harmonise_data_female)
mr_singlesnp_combined <- mr_singlesnp(harmonise_data_combined)

mr_hetrogeneity_male <- mr_heterogeneity(harmonise_data_male, method_list = methods_heterogeneity)
mr_hetrogeneity_female <- mr_heterogeneity(harmonise_data_female, method_list = methods_heterogeneity)
mr_hetrogeneity_combined <- mr_heterogeneity(harmonise_data_combined, method_list = methods_heterogeneity)

mr_pleiotropy_male <- mr_pleiotropy_test(harmonise_data_male)
mr_pleiotropy_female <- mr_pleiotropy_test(harmonise_data_female)
mr_pleiotropy_combined <- mr_pleiotropy_test(harmonise_data_combined)

mr_leaveoneout_male <- mr_leaveoneout(harmonise_data_male)
mr_leaveoneout_female <- mr_leaveoneout(harmonise_data_female)
mr_leaveoneout_combined <- mr_leaveoneout(harmonise_data_combined)

## Plots ====
source("scripts/my_mr_scatter_plot.R")
plot_mr_scatter_male <- my_mr_scatter_plot(mr_results_male, harmonise_data_male)
plot_mr_scatter_female <- my_mr_scatter_plot(mr_results_female, harmonise_data_female)
plot_mr_scatter_combined <- my_mr_scatter_plot(mr_results_combined, harmonise_data_combined)

plot_singlesnp_forest_male <- mr_forest_plot(mr_singlesnp_male)
plot_singlesnp_forest_female <- mr_forest_plot(mr_singlesnp_female)
plot_singlesnp_forest_combined <- mr_forest_plot(mr_singlesnp_combined)

plot_leaveoneout_forest_male <- mr_leaveoneout_plot(mr_leaveoneout_male)
plot_leaveoneout_forest_female <- mr_leaveoneout_plot(mr_leaveoneout_female)
plot_leaveoneout_forest_combined <- mr_leaveoneout_plot(mr_leaveoneout_combined)

plot_mr_funnel_male <- mr_funnel_plot(mr_singlesnp_male)
plot_mr_funnel_female <- mr_funnel_plot(mr_singlesnp_female)
plot_mr_funnel_combined <- mr_funnel_plot(mr_singlesnp_combined)

### save plots ====
pdf("analysis/001_adiposity_colorectal/figures/plot_mr_scatter_male.pdf")
for (i in 1:length(plot_mr_scatter_male)) {
  print(plot_mr_scatter_male[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_mr_scatter_female.pdf")
for (i in 1:length(plot_mr_scatter_female)) {
  print(plot_mr_scatter_female[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_mr_scatter_combined.pdf")
for (i in 1:length(plot_mr_scatter_combined)) {
  print(plot_mr_scatter_combined[[i]])
}
dev.off()

pdf("analysis/001_adiposity_colorectal/figures/plot_singlesnp_forest_male.pdf")
for (i in 1:length(plot_singlesnp_forest_male)) {
  print(plot_singlesnp_forest_male[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_singlesnp_forest_female.pdf")
for (i in 1:length(plot_singlesnp_forest_female)) {
  print(plot_singlesnp_forest_female[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_singlesnp_forest_combined.pdf")
for (i in 1:length(plot_singlesnp_forest_combined)) {
  print(plot_singlesnp_forest_combined[[i]])
}
dev.off()

pdf("analysis/001_adiposity_colorectal/figures/plot_leaveoneout_forest_male.pdf")
for (i in 1:length(plot_leaveoneout_forest_male)) {
  print(plot_leaveoneout_forest_male[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_leaveoneout_forest_female.pdf")
for (i in 1:length(plot_leaveoneout_forest_female)) {
  print(plot_leaveoneout_forest_female[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_leaveoneout_forest_combined.pdf")
for (i in 1:length(plot_leaveoneout_forest_combined)) {
  print(plot_leaveoneout_forest_combined[[i]])
}
dev.off()

pdf("analysis/001_adiposity_colorectal/figures/plot_mr_funnel_male.pdf")
for (i in 1:length(plot_mr_funnel_male)) {
  print(plot_mr_funnel_male[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_mr_funnel_female.pdf")
for (i in 1:length(plot_mr_funnel_female)) {
  print(plot_mr_funnel_female[[i]])
}
dev.off()
pdf("analysis/001_adiposity_colorectal/figures/plot_mr_funnel_combined.pdf")
for (i in 1:length(plot_mr_funnel_combined)) {
  print(plot_mr_funnel_combined[[i]])
}
dev.off()

## Save output ====
write.table(exposure_data, "analysis/001_adiposity_colorectal/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(outcome_data_combined, "analysis/001_adiposity_colorectal/outcome_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(harmonise_data_male, "analysis/001_adiposity_colorectal/harmonise_data_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(harmonise_data_female, "analysis/001_adiposity_colorectal/harmonise_data_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(harmonise_data_combined, "analysis/001_adiposity_colorectal/harmonise_data_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_results_male, "analysis/001_adiposity_colorectal/mr_results_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_results_female, "analysis/001_adiposity_colorectal/mr_results_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_results_combined, "analysis/001_adiposity_colorectal/mr_results_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_singlesnp_male, "analysis/001_adiposity_colorectal/mr_singlesnp_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_singlesnp_female, "analysis/001_adiposity_colorectal/mr_singlesnp_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_singlesnp_combined, "analysis/001_adiposity_colorectal/mr_singlesnp_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_hetrogeneity_male, "analysis/001_adiposity_colorectal/mr_hetrogeneity_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_hetrogeneity_female, "analysis/001_adiposity_colorectal/mr_hetrogeneity_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_hetrogeneity_combined, "analysis/001_adiposity_colorectal/mr_hetrogeneity_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_pleiotropy_male, "analysis/001_adiposity_colorectal/mr_pleiotropy_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_pleiotropy_female, "analysis/001_adiposity_colorectal/mr_pleiotropy_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_pleiotropy_combined, "analysis/001_adiposity_colorectal/mr_pleiotropy_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_leaveoneout_male, "analysis/001_adiposity_colorectal/mr_leaveoneout_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_leaveoneout_female, "analysis/001_adiposity_colorectal/mr_leaveoneout_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(mr_leaveoneout_combined, "analysis/001_adiposity_colorectal/mr_leaveoneout_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
