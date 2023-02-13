rm(list=ls())

# environment ====
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#remotes::install_github("MRCIEU/genetics.binaRies", force = F)
#remotes::install_github("explodecomputer/plinkbinr", force = F)
#remotes::install_github("chr1swallace/coloc@main", force = F)
#remotes::install_github("sjmgarnier/viridis", force = F)
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)

source("scripts/011_colocalisation/functions/my_coloc_chriswallace.R")

# data ====
a <- read.table("analysis/008_mvmr/mvmr_results.txt", header = T, sep = "\t")
a <- subset(a, group == "Female")
a <- subset(a, exposure != "BMI")
a <- subset(a, exposure != "WHR")
a <- subset(a, exposure != "WHRadjBMI")
a <- subset(a, cancer == "rectal")

filenames_mvmr <- paste0("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb//", a$exposure, "_", a$gene, "_", a$protein, ".txt.gz.annotated.gz.exclusions.gz.alleles.gz.unzipped.cis.txt")
filenames_all <- dir("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb/", recursive = TRUE, full.names = TRUE, pattern = ".cis.txt")
filenames <- intersect(filenames_mvmr, filenames_all)

# exposure ====
exposure_list <- lapply(filenames, fread, col.names = c("CHR", "POS", "SNPID", "SNP", "EA", "OA",
                                                        "beta.exposure", "pval.exposure", "minus_log10_pval", "se.exposure", "samplesize.exposure",
                                                        "EAF", "exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure"))
length(exposure_list)
exposure_list <- exposure_list[sapply(exposure_list, nrow) > 0]
length(exposure_list)
exposure_list <- purrr::discard(exposure_list, ~any(.x$CHR == "chrX")) # X CHR not available in outcome
length(exposure_list)

# format exposure ====
exposure_filenames <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb//", "", filenames)
exposure_filenames <- gsub(".txt.gz.annotated.gz.exclusions.gz.alleles.gz.unzipped.cis.txt", "", exposure_filenames)

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$exposure <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/", "", exposure_list[[i]]$exposure)
}

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$exposure <- gsub(".txt.gz.unzipped", "", exposure_list[[i]]$exposure)
}

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$id.exposure <- paste0(exposure_filenames[[i]], "_", "joint_rectal_female")
}

# outcome data ====
filenames <- c("joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt")

outcome_list <- list()

for (i in 1:length(exposure_list)){
  
  outcome_list[i] <- lapply(paste0("/data/GWAS_data/files/huyghe_2018_PMID30510241/processed/",filenames), 
                            read_outcome_data,
                            snps = exposure_list[[i]]$SNP,
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
                            phenotype_col = "phenotype")
  
  outcome_list[[i]]$outcome <- "joint_rectal_female"
  outcome_list[[i]]$id.outcome <- paste0(exposure_filenames[[i]], "_", outcome_list[[i]]$outcome)
  
}

# harmonise ====
exposure <- bind_rows(exposure_list)
outcome <- bind_rows(outcome_list)
harmonise_data <- harmonise_data(exposure, outcome, action = 2)
harmonise_data$remove_duplicates <- paste0(harmonise_data$SNP, "_", harmonise_data$id.exposure)
harmonise_data <- harmonise_data[!duplicated(harmonise_data$remove_duplicates),]
harmonise_data_list <- split(harmonise_data, harmonise_data$id.exposure)

# loop over all harmonised data and run ld matrix, formatting, coloc, save ====
table_master <- data.frame() # make empty dataframe for final results

for (i in 1:length(harmonise_data_list)){
  
  label_exposure <- unique(harmonise_data_list[[i]]$exposure)
  label <- paste0(label_exposure, "_", "joint_rectal_female")
  label_outcome <- "joint_rectal_female"
  
  # make ld matrix ====
  ld <- ld_matrix_local(
    harmonise_data_list[[i]]$SNP,
    with_alleles = FALSE, 
    bfile = "/data/GWAS_data/files/references/1kG_v3/EUR/EUR",
    plink_bin = get_plink_exe())
  
  # format LD matrix and harmonised list ====
  ld <- ld[which(rownames(ld) %in% harmonise_data_list[[i]]$SNP), which(colnames(ld) %in% harmonise_data_list[[i]]$SNP)]
  harmonise_data_list[[i]] <- harmonise_data_list[[i]][which(harmonise_data_list[[i]]$SNP %in% rownames(ld)),]
  ld <- ld[match(harmonise_data_list[[i]]$SNP,rownames(ld)),]
  ld <- ld[,match(harmonise_data_list[[i]]$SNP, colnames(ld))]
  harmonise_data_list[[i]] <- harmonise_data_list[[i]][match(rownames(ld), harmonise_data_list[[i]]$SNP),]
  
  # make lists for coloc ====
  coloc_data_exposure <- list(beta = harmonise_data_list[[i]]$beta.exposure, varbeta = harmonise_data_list[[i]]$se.exposure^2, MAF = harmonise_data_list[[i]]$eaf.exposure, type = "quant", N = 35559, snp = rownames(ld), LD = ld, position = harmonise_data_list[[i]]$POS)
  coloc_data_outcome <- list(beta = harmonise_data_list[[i]]$beta.outcome, varbeta = harmonise_data_list[[i]]$se.outcome^2, MAF = harmonise_data_list[[i]]$eaf.outcome, type = "cc", N = 120328, snp = rownames(ld), LD = ld, position = harmonise_data_list[[i]]$POS)
  
  # coloc ====  
  coloc_results <- coloc.abf(dataset1 = coloc_data_exposure, dataset2 = coloc_data_outcome)
  
  pdf(paste0("analysis/009_colocalisation/results/joint_rectal_female/figures/", label, ".pdf"), 
      height = 10, width = 10)
  coloc_sensitivity <- my_sensitivity(coloc_results, "H4 > 0.9", 
                                      trait1_title = label_exposure, trait2_title = label_outcome)
  dev.off()
  
  # save ====
  saveRDS(coloc_results, paste0("analysis/009_colocalisation/results/joint_rectal_female/", label, ".RData"))
  
  # make table ====
  table <- data.frame(
    exposure = label_exposure,
    outcome = label_outcome,
    id = label,
    nsnps = coloc_results["summary"][[1]][1],
    h0 = coloc_results["summary"][[1]][2],
    h1 = coloc_results["summary"][[1]][3],
    h2 = coloc_results["summary"][[1]][4],
    h3 = coloc_results["summary"][[1]][5],
    h4 = coloc_results["summary"][[1]][6])
  
  table_master <- rbind(table_master, table)
  
}

write.table(table_master, "analysis/009_colocalisation/results/joint_rectal_female/001_coloc_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
