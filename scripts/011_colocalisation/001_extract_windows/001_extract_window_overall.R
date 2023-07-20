rm(list=ls())

# environment ====
library(data.table)
library(TwoSampleMR)
library(dplyr)

## cis/window data ====
filenames <- dir("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb/", recursive = TRUE, full.names = TRUE, pattern = ".cis.txt")

exposure_list <- lapply(filenames, fread, col.names = c("CHR", "POS", "SNPID", "SNP", "EA", "OA",
                                                        "beta.exposure", "pval.exposure", "minus_log10_pval", "se.exposure", "samplesize.exposure",
                                                        "EAF", "exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure"))

length(exposure_list)
exposure_list <- exposure_list[sapply(exposure_list, nrow) > 0]
length(exposure_list)
exposure_list <- purrr::discard(exposure_list, ~any(.x$CHR == "chrX")) # X CHR not available in outcome
length(exposure_list)

exposure_filenames <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb//", "", filenames)
exposure_filenames <- gsub(".txt.gz.annotated.gz.exclusions.gz.alleles.gz.unzipped.cis.txt", "", exposure_filenames)

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$exposure <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/", "", exposure_list[[i]]$exposure)
}

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$exposure <- gsub(".txt.gz.unzipped", "", exposure_list[[i]]$exposure)
}

for (i in 1:length(exposure_list)){
  exposure_list[[i]]$id.exposure <- paste0(exposure_filenames[[i]], "_", "overall")
}

## extract outcome data ====
filenames <- c("MarginalMeta_125K_Results.tsv.annotated.txt")

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
    
    outcome_list[[i]]$outcome <- "overall"
    outcome_list[[i]]$id.outcome <- paste0(exposure_filenames[[i]], "_", outcome_list[[i]]$outcome)

}

outcomes <- bind_rows(outcome_list)

write.table(outcomes, "analysis/009_colocalisation/outcome_data/000_overall_new.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
