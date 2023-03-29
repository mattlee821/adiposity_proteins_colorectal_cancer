rm(list=ls())
## set environment ====
library(dplyr)

# adiposity ====
filenames <- dir("data/000_mvmr_data/adiposity/", recursive = TRUE, full.names = TRUE, pattern = ".txt")
list <- lapply(filenames, read.table, header = F, sep = " ", fill = T, col.names = c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele", "BETA", "SE", "P", "N", "INFO",	"phenotype")) 
data <- bind_rows(list)
data <- data %>%
  separate(INFO, c("INFO", "phenotype"), "\t")

# sort names out 
table(data$phenotype)
data$group[data$phenotype == "bmi_combined.txt"] <- "Sex-combined"
data$group[data$phenotype == "whr_combined.txt"] <- "Sex-combined"

data$group[data$phenotype == "bmi_male.txt"] <- "Male"
data$group[data$phenotype == "whr_male.txt"] <- "Male"

data$group[data$phenotype == "bmi_female.txt"] <- "Female"
data$group[data$phenotype == "whr_female.txt"] <- "Female"

data$adiposity[data$phenotype == "bmi_combined.txt"] <- "BMI"
data$adiposity[data$phenotype == "whr_combined.txt"] <- "WHR"

data$adiposity[data$phenotype == "bmi_male.txt"] <- "BMI"
data$adiposity[data$phenotype == "whr_male.txt"] <- "WHR"

data$adiposity[data$phenotype == "bmi_female.txt"] <- "BMI"
data$adiposity[data$phenotype == "whr_female.txt"] <- "WHR"

instruments <- data

# associations ====
associations <- read.table("analysis/008_mvmr/000_snp_list.txt", header = T, sep = "\t")

# join ====
data <- left_join(instruments, associations, by = c("adiposity", "group", "SNP"))
data <- data[complete.cases(data), ]

colnames(data) <- c(
  "CHR", "POS", "SNP", "effect_allele.exposure", "other_allele.exposure",
  "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", 
  "INFO", "phenotype", "group", "adiposity", 
  "id.exposure", "SNP_from", "sequence_id", "cancer")

data <- data[,c("CHR", "POS", "SNP", "effect_allele.exposure", "other_allele.exposure",
                "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                "id.exposure", "adiposity","sequence_id", "cancer", "SNP_from", "group")]

data$exposure <- data$id.exposure 

# save
write.table(data, "data/000_mvmr_data/001_adiposity_data.txt", 
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# proteins ====
filenames <- dir("data/000_mvmr_data/proteins/", recursive = TRUE, full.names = TRUE, pattern = ".txt")
list <- lapply(filenames, read.table, header = F, sep = "\t", fill = T, col.names = c("Chrom", "Pos", "Name", "SNP", "effectAllele", "otherAllele", "Beta", "pval", "min_log10_pval", "SE", "N", "ImpMAF", "phenotype", "EA", "OA", "EAF")) 
data <- bind_rows(list)

# format protein names
table(data$phenotype)
data$phenotype <- gsub(".txt.gz.unzipped", "", data$phenotype)
data$phenotype <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/processed/", "", data$phenotype)

# sequence_id column 
ncols <- max(stringr::str_count(data$phenotype, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = phenotype,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

# protein_name column 
data$phenotype <- sub(".*?_", "", data$phenotype) # remove first part of sequence_id from column
data$phenotype <- sub(".*?_", "", data$phenotype) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$phenotype) # remove genename

# gene name column 
data$gene <- gsub("\\_.*","",data$phenotype)

# remove extra cols 
instruments <- data[, -which(names(data) %in% c("col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"))]

# associations ====
associations <- read.table("analysis/008_mvmr/000_snp_list.txt", header = T, sep = "\t")

# join ====
data <- left_join(instruments, associations, by = c("sequence_id", "SNP"))
data <- data[complete.cases(data), ]

data <- data[,c("Chrom", "Pos", "SNP", "EA", "OA",
                "EAF", "Beta", "SE", "pval", "N",
                "id.exposure", "adiposity","sequence_id", "cancer", "exposure", "group")]

colnames(data) <- c(
  "CHR", "POS", "SNP", "effect_allele.outcome", "other_allele.outcome",
  "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", 
  "id.outcome", "adiposity", "sequence_id", "cancer", "SNP_from", "group")

data$outcome <- data$id.outcome

# save
write.table(data, "data/000_mvmr_data/001_protein_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# cancer ====
filenames <- dir("data/000_mvmr_data/cancer/", recursive = TRUE, full.names = TRUE, pattern = ".txt")
filenames <- filenames[-10]
list <- lapply(filenames, read.table, header = F)
list1 <- list[c(1:8)]
list2 <- list[c(9:15)]

data1 <- bind_rows(list1)
data1 <- data1[,c(17,18,19,2,3,4,8,9,10,25)]
colnames(data1) <- c("CHR","POS", "SNP", "effectAllele.outcome","otherAllele.outcome","eaf.outcome","beta.outcome","se.outcome","pval.outcome","outcome")
data2 <- bind_rows(list2)
data2 <- data2[,c(2,3,25,4,5,6,10,11,12,31)]
colnames(data2) <- c("CHR","POS", "SNP", "effectAllele.outcome","otherAllele.outcome","eaf.outcome","beta.outcome","se.outcome","pval.outcome","outcome")

data <- bind_rows(data1,data2)

# format  names
table(data$outcome)
data$group[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "Female"
data$group[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "Male"
data$group[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- "Sex-combined"
data$group[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- "Female"
data$group[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- "Male"

data$cancer[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$cancer[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- "colon"
data$cancer[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$cancer[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "distal"
data$cancer[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$cancer[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "proximal"
data$cancer[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$cancer[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- "rectal"
data$cancer[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- "overall"
data$cancer[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- "colon"
data$cancer[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- "distal"
data$cancer[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- "proximal"
data$cancer[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- "rectal"
data$cancer[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- "overall"
data$cancer[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- "overall"

# add N for each CRC outcome
data$samplesize.outcome[data$outcome == "MarginalMeta_125K_Results.tsv.annotated.txt"] <- (58131 + 67347)

data$samplesize.outcome[data$outcome == "Stratified_colon_Meta_125K_Results.tsv.annotated.txt"] <- (32002 + 64159)
data$samplesize.outcome[data$outcome == "Stratified_distal_Meta_125K_Results.tsv.annotated.txt"] <- (14376 + 64159)
data$samplesize.outcome[data$outcome == "Stratified_proximal_Meta_125K_Results.tsv.annotated.txt"] <- (15706 + 64159)
data$samplesize.outcome[data$outcome == "Stratified_rectal_Meta_125K_Results.tsv.annotated.txt"] <- (16212 + 64159)

data$samplesize.outcome[data$outcome == "Stratified_male_Meta_125K_Results.tsv.annotated.txt"] <- (31288 + 34527)
data$samplesize.outcome[data$outcome == "joint_colon_Male_wald_MAC50_1.TBL.annotated.txt"] <- 34527
data$samplesize.outcome[data$outcome == "joint_distal_Male_wald_MAC50_1.TBL.annotated.txt"] <- 34527
data$samplesize.outcome[data$outcome == "joint_proximal_Male_wald_MAC50_1.TBL.annotated.txt"] <- 34527
data$samplesize.outcome[data$outcome == "joint_rectal_Male_wald_MAC50_1.TBL.annotated.txt"] <- 34527

data$samplesize.outcome[data$outcome == "Stratified_female_Meta_125K_Results.tsv.annotated.txt"] <- (26834 + 32820)
data$samplesize.outcome[data$outcome == "joint_colon_Female_wald_MAC50_1.TBL.annotated.txt"] <- 32820
data$samplesize.outcome[data$outcome == "joint_distal_Female_wald_MAC50_1.TBL.annotated.txt"] <- 32820
data$samplesize.outcome[data$outcome == "joint_proximal_Female_wald_MAC50_1.TBL.annotated.txt"] <- 32820
data$samplesize.outcome[data$outcome == "joint_rectal_Female_wald_MAC50_1.TBL.annotated.txt"] <- 32820

instruments <- data

# associations ====
associations <- read.table("analysis/008_mvmr/000_snp_list.txt", header = T, sep = "\t")

# join ====
data <- left_join(instruments, associations, by = c("cancer", "group", "SNP"))
data <- data[complete.cases(data), ]

data <- data[,c("CHR", "POS", "SNP", "effectAllele.outcome", "otherAllele.outcome",
                "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", 
                "id.exposure", "adiposity","sequence_id", "cancer", "cancer", "group")]

colnames(data) <- c(
  "CHR", "POS", "SNP", "effect_allele.outcome", "other_allele.outcome",
  "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", 
  "id.outcome", "adiposity", "sequence_id", "cancer", "SNP_from", "group")

data$outcome <- data$id.outcome

# save
write.table(data, "data/000_mvmr_data/001_cancer_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
