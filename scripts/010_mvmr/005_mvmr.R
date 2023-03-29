rm(list=ls())
# MR analysis of measures of adiposity and metabolites 

# environment ====
## library ====
library(TwoSampleMR)
library(MVMR)
library(data.table)
library(dplyr)

# data ====
adiposity_proteins <- read.table("analysis/008_mvmr/harmonise_data_adiposity_proteins.txt", header = T, sep = "\t")
adiposity_cancer <- read.table("analysis/008_mvmr/harmonise_data_adiposity_cancer.txt", header = T, sep = "\t")

# format ====
adiposity_proteins <- adiposity_proteins[,c(
    "SNP", 
    "effect_allele.exposure", "other_allele.exposure",
    "effect_allele.outcome", "other_allele.outcome",
    "beta.exposure","beta.outcome", 
    "eaf.exposure", "eaf.outcome",
    "se.exposure","se.outcome",
    "pval.exposure","pval.outcome",
    "samplesize.exposure","samplesize.outcome",
    "exposure","id.exposure",
    "id.outcome","outcome")]

colnames(adiposity_proteins) <- c(
    "SNP", 
    "effect_allele.exposure", "other_allele.exposure",
    "effect_allele.outcome1", "other_allele.outcome1",
    "beta.exposure","beta.outcome1", 
    "eaf.exposure", "eaf.outcome1",
    "se.exposure","se.outcome1",
    "pval.exposure","pval.outcome1",
    "samplesize.exposure","samplesize.outcome1",
    "exposure","id.exposure",
    "id.outcome1","outcome1")

adiposity_cancer <- adiposity_cancer[,c(
  "SNP", 
  "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome",
  "beta.exposure","beta.outcome", 
  "eaf.exposure", "eaf.outcome",
  "se.exposure","se.outcome",
  "pval.exposure","pval.outcome",
  "samplesize.exposure","samplesize.outcome",
  "exposure","id.exposure",
  "id.outcome","outcome")]

colnames(adiposity_cancer) <- c(
  "SNP", 
  "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome1", "other_allele.outcome1",
  "beta.exposure","beta.outcome", 
  "eaf.exposure", "eaf.outcome",
  "se.exposure","se.outcome",
  "pval.exposure","pval.outcome",
  "samplesize.exposure","samplesize.outcome",
  "exposure","id.exposure",
  "id.outcome","outcome")

#make harmonised data frame ====
harmonise_data <- left_join(adiposity_proteins, adiposity_cancer, by = c("SNP", "id.exposure"))
harmonise_data <- harmonise_data[,c("SNP",
                                    "beta.exposure.x","beta.outcome1", "beta.outcome",
                                    "se.exposure.x","se.outcome1","se.outcome","outcome")]

colnames(harmonise_data) <- c(
    "SNP",
    "adiposity_b", "protein_b", "cancer_b",
    "adiposity_se","protein_se", "cancer_se",
    "ID")

harmonise_data <- subset(harmonise_data, !grepl("HRC-EUR", harmonise_data[,"ID"]))

harmonise_outcome <- split(harmonise_data, f = harmonise_data$ID) # make a list of data frames with each data frame for each ID

# make MVMR data ====
mvmr_data <- list()
for (i in 1:length(harmonise_outcome))
  mvmr_data[[i]] <- format_mvmr(
    BXGs = harmonise_outcome[[i]][,c(2,3)],
    BYG = harmonise_outcome[[i]]$cancer_b,
    seBXGs = harmonise_outcome[[i]][,c(5,6)],
    seBYG = harmonise_outcome[[i]]$cancer_se, 
    RSID = harmonise_outcome[[i]]$SNP)

# run MVMR ====
result <- list()
for (i in 1:length(mvmr_data))
  result[[i]] <- as.data.frame(ivw_mvmr(r_input = mvmr_data[[i]]))

## run sensitivity
fstat <- list()
for (i in 1:length(mvmr_data))
  fstat[[i]] <- strength_mvmr(r_input = mvmr_data[[i]], gencov = 0)

qstat <- list()
for (i in 1:length(mvmr_data))
  qstat[[i]] <- pleiotropy_mvmr(r_input = mvmr_data[[i]], gencov = 0)

# make results table ====
table <- list()
for (i in 1:length(result))
  table[[i]] <- data.frame(
    b = c(result[[i]][1,1], result[[i]][2,1]),
    se = c(result[[i]][1,2], result[[i]][2,2]),
    t = c(result[[i]][1,3], result[[i]][2,3]),
    p = c(result[[i]][1,4], result[[i]][2,4]),
    ID = c(levels(as.factor(harmonise_data$ID))[[i]], levels(as.factor(harmonise_data$ID))[[i]]),
    fstat = c(fstat[[i]]),
    qstat = c(qstat[[i]][1]),
    qstat_p = c(qstat[[i]][2]))

# format results table ====
## sequence_id column 
ncols <- max(stringr::str_count(table[[1]]$ID, "_")) + 1
colmn <- paste0("col", 1:ncols)
for (i in 1:length(table))
  table[[i]] <- 
  tidyr::separate(
    data = table[[i]],
    col = ID,
    sep = "_",
    into = colmn,
    remove = FALSE)

for (i in 1:length(table))
  table[[i]]$adiposity <- table[[i]]$col1

for (i in 1:length(table))
  table[[i]]$sequence_id <- paste0(table[[i]]$col2, "_", table[[i]]$col3)

for (i in 1:length(table))
  table[[i]]$cancer <- table[[i]]$col4

for (i in 1:length(table))
  table[[i]]$group <- table[[i]]$col5

for (i in 1:length(table))
  table[[i]] <- table[[i]][, -which(names(table[[i]]) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"))]

mvmr_results <- list()
for (i in 1:length(table))
  mvmr_results[[i]] <- data.frame(
    b = table[[i]]$b,
    se = table[[i]]$se,
    t = table[[i]]$t,
    p = table[[i]]$p,
    ID = table[[i]]$ID,
    fstat_adiposity = table[[i]]$fstat.exposure1,
    fstat_protein = table[[i]]$fstat.exposure2,
    Qstat = table[[i]]$Qstat,
    Qpval = table[[i]]$Qpval,
    adiposity = table[[i]]$adiposity,
    sequence_ID = table[[i]]$sequence_id,
    cancer = table[[i]]$cancer,
    group = table[[i]]$group,
    exposure = c(table[[i]][1,10], table[[i]][2,11]),
    mediator = c(table[[i]][1,11], table[[i]][2,10])
  )

mvmr_results <- bind_rows(mvmr_results)

# OR and CI
mvmr_results$OR <- exp(mvmr_results$b)
mvmr_results$lower_ci <- exp(mvmr_results$b - (1.96 * mvmr_results$se))
mvmr_results$upper_ci <- exp(mvmr_results$b + (1.96 * mvmr_results$se))

# join protein info ====
proteins <- readxl::read_xlsx("data/protein_info.xlsx")
mvmr_results <- left_join(mvmr_results, proteins, by = "sequence_ID")

# rearrange ====
mvmr_results <- mvmr_results[,c(
  "ID", "b", "se", "OR", "lower_ci", "upper_ci", "p", "Qstat", "Qpval",
  "fstat_adiposity", "fstat_protein",
  "exposure", "mediator", "cancer", "group", "adiposity", "protein", 
  "protein_name", "gene", "uni_prot", "ensembl_gene_ID", "entrez_gene_ID", "HGNC_ID"
)]

# save ====
write.table(mvmr_results, "analysis/008_mvmr/mvmr_results.txt",
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


