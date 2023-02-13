rm(list=ls())

# source ====

# data ====
data <- read.table("analysis/008_mvmr/mvmr_results.txt", header = T, sep = "\t")
labels <- unique(data$ID)
adiposity_labels <- c("BMI", "WHR")

adiposity_mvmr <- data[data$exposure %in% adiposity_labels, ]

# adiposity_cancer ====
data <- read.table("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T, sep = "\t")
colnames(data) <- c("exposure","cancer", "group", "nsnp", "method", "b_MR", "se_MR", "pval_MR", "OR_MR", "lower_ci_MR", "upper_ci_MR")
adiposity_cancer <- subset(data, method == "IVW-MRE")

# join data ====
data <- left_join(adiposity_mvmr, adiposity_cancer, by = c("exposure", "cancer", "group"))

# proportion mediated
data$total_effect <- data$b_MR
data$direct_effect <- data$b
data$indirect_effect <- data$total_effect - data$direct_effect
data$proportion_mediated <- (data$indirect_effect / data$total_effect) * 100

data <- data[,c(1,32:35)]

write.table(data, "analysis/008_mvmr/proportion_mediated.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")






# univariable MR proportion mediated ====
# adiposity_cancer ====
data <- read.table("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T, sep = "\t")
colnames(data) <- c("adiposity","cancer", "group", "nsnp", "method", "b_adiposity_cancer_MR", "se_adiposity_cancer_MR", "pval_MR", "OR_MR", "lower_ci_MR", "upper_ci_MR")
data <- data[,c("adiposity","cancer", "group","method","b_adiposity_cancer_MR")]
adiposity_cancer <- subset(data, method == "IVW-MRE")

# adiposity_protein ====
data <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")
colnames(data) <- c("adiposity","group", "outcome", "sequence_id", "protein", "gene", "nsnp", "method", "b_adiposity_protein_MR", "se_adiposity_protein_MR", "pval_MR", "lower_ci_MR", "upper_ci_MR")
data <- data[,c("adiposity","group","outcome", "sequence_id", "protein", "gene","method","b_adiposity_protein_MR")]
adiposity_protein <- subset(data, method == "IVW-MRE")

# protein_cancer ====
data <- read.table("analysis/003_proteins_colorectal/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
colnames(data) <- c("exposure", "group", "cancer", "sequence_id", "protein", "gene", "nsnp", "method", "b_protein_cancer_MR", "se_protein_cancer_MR", "pval_MR", "lower_ci_MR", "upper_ci_MR")
data <- data[,c("exposure","group", "sequence_id", "protein", "gene", "cancer", "method", "b_protein_cancer_MR")]
protein_cancer <- subset(data, method == "IVW-MRE")

# associations ====
data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/004_adiposity_protein_cancer_pos_neg.txt", header = T, sep = "\t")
names(data)
associations <- data[,c("adiposity","sequence_id","protein","gene","cancer","group")]

data <- left_join(associations, adiposity_cancer, by = c("adiposity","cancer", "group"))
names(data)

data1 <- left_join(data, adiposity_protein, by = c("adiposity","protein","group"))
names(data1)
data1 <- data1[,c("adiposity","sequence_id.x","protein","gene.x","cancer","group", "method.x","b_adiposity_cancer_MR", "b_adiposity_protein_MR")]
colnames(data1) <- c("adiposity","sequence_id","protein","gene","cancer","group", "method","b_adiposity_cancer_MR", "b_adiposity_protein_MR")
data1 <- unique(data1)

data1 <- left_join(data1, protein_cancer, by = c("group", "sequence_id", "protein", "gene", "cancer"))
data1 <- data1[,c("adiposity","sequence_id","protein","gene","cancer","group", "method.x","b_adiposity_cancer_MR", "b_adiposity_protein_MR","b_protein_cancer_MR")]
colnames(data1) <- c("adiposity","sequence_id","protein","gene","cancer","group", "method","b_adiposity_cancer_MR", "b_adiposity_protein_MR","b_protein_cancer_MR")

data1$ID <- paste0(data1$adiposity, "_", data1$sequence_id, "_", data1$cancer, "_", data1$group)

a <- unique(data1$ID)

data1$total_effect <- data1$b_adiposity_cancer_MR
data1$indirect_effect <- data1$b_adiposity_protein_MR * data1$b_protein_cancer_MR
data1$direct_effect <- data1$total_effect - data1$indirect_effect

data1$proportion_mediated <- (data1$indirect_effect / data1$total_effect) * 100
a <- subset(data1, duplicated(ID))



data <- left_join(adiposity_mvmr, adiposity_cancer, by = c("exposure", "cancer", "group"))

