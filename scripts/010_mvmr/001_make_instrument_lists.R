rm(list=ls())
# MR analysis of measures of adiposity and metabolites 

# environment ====
## library ====
library(TwoSampleMR)
library(data.table)
library(dplyr)

# associated list ====
associations <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/004_adiposity_protein_cancer_pos_neg.txt", header = T, sep = "\t")
associated_proteins <- unique(associations[,c("sequence_id")])
associations <- associations[,c("adiposity", "sequence_id", "cancer", "group")]

# adiposity ====
data <- read.table("data/exposure_data/000_clumped_adiposity.txt", header = T, sep = "\t")

# sort names out 
table(data$id.exposure)
data$group[data$id.exposure == "BMI_combined"] <- "Sex-combined"
data$group[data$id.exposure == "WHR_combined"] <- "Sex-combined"

data$group[data$id.exposure == "BMI_male"] <- "Male"
data$group[data$id.exposure == "WHR_male"] <- "Male"

data$group[data$id.exposure == "BMI_female"] <- "Female"
data$group[data$id.exposure == "WHR_female"] <- "Female"


data$adiposity[data$id.exposure == "BMI_combined"] <- "BMI"
data$adiposity[data$id.exposure == "WHR_combined"] <- "WHR"

data$adiposity[data$id.exposure == "BMI_male"] <- "BMI"
data$adiposity[data$id.exposure == "WHR_male"] <- "WHR"

data$adiposity[data$id.exposure == "BMI_female"] <- "BMI"
data$adiposity[data$id.exposure == "WHR_female"] <- "WHR"

adiposity <- data[,c("exposure", "adiposity", "group", "SNP")]

# proteins ====
data <- read.table("data/000_protein_data/exposure_data/000_clumped.txt", header = T, sep = "\t")
proteins <- data[,c("exposure", "sequence_id", "SNP")]
proteins <- proteins[proteins$sequence_id %in% associated_proteins, ]

# make lists ====
a <- left_join(adiposity, associations, by = c("adiposity", "group"))
a <- a[complete.cases(a), ]
b <- left_join(proteins, associations, by = "sequence_id")
b <- b[complete.cases(b), ]

data <- rbind(a,b)
data$id.exposure <- paste0(data$adiposity, "_", data$sequence_id, "_", data$cancer, "_", data$group)
data <- data[,c("id.exposure", "SNP", "exposure", "adiposity", "sequence_id", "cancer", "group")]

data_clump <- clump_data(data,
           clump_kb = 10000,
           clump_r2 = 0.001)

write.table(data, "analysis/008_mvmr/000_snp_list.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

a <- as.data.frame(unique(data$SNP))

write.table(a, "analysis/008_mvmr/000_unique_snp_list.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
