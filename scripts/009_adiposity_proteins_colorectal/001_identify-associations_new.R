rm(list=ls())
library(dplyr)

# adiposity-cancer associations ====
## adiposity-cancer
adiposity_cancer <- read.table("analysis/001_adiposity_colorectal/001_directional_consistency.txt", header = T, sep = "\t")
adiposity_cancer <- subset(adiposity_cancer, outcome != "overall-HRC-EUR")
table(adiposity_cancer$outcome, adiposity_cancer$group, adiposity_cancer$exposure)
adiposity_cancer$id <- paste0(adiposity_cancer$exposure, "_", adiposity_cancer$outcome, "_", adiposity_cancer$group, "_", adiposity_cancer$direction_group)
adiposity_cancer_id <- adiposity_cancer$id

## adiposity-cancer
cancer_adiposity <- read.table("analysis/004_colorectal_adiposity/001_directional_consistency.txt", header = T, sep = "\t")
cancer_adiposity <- subset(cancer_adiposity, exposure != "overall-HRC-EUR")
cancer_adiposity$id <- paste0(cancer_adiposity$outcome, "_", cancer_adiposity$exposure, "_", cancer_adiposity$group, "_", cancer_adiposity$direction_group)
cancer_adiposity_id <- cancer_adiposity$id

## remove inconsistencies
a <- intersect(adiposity_cancer_id, cancer_adiposity_id)
adiposity_cancer <- adiposity_cancer[!adiposity_cancer$id %in% a, ]
table(adiposity_cancer$outcome, adiposity_cancer$group, adiposity_cancer$exposure)

# adiposity-protein associations ====
## adiposity-protein
adiposity_proteins <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
adiposity_proteins <- adiposity_proteins[,c("sequence_id", "protein", "gene", "exposure", "group", "b")]
colnames(adiposity_proteins)[4] <- "adiposity"

table1 <- data.frame(
  analysis = c("adiposity > protein"),
  
  bmi_combined = c(nrow(subset(adiposity_proteins, adiposity == "BMI" & group == "Sex-combined"))),
  bmi_male = c(nrow(subset(adiposity_proteins, adiposity == "BMI" & group == "Male"))),
  bmi_female = c(nrow(subset(adiposity_proteins, adiposity == "BMI" & group == "Female"))),
  
  whr_combined = c(nrow(subset(adiposity_proteins, adiposity == "WHR" & group == "Sex-combined"))),
  whr_male = c(nrow(subset(adiposity_proteins, adiposity == "WHR" & group == "Male"))),
  whr_female = c(nrow(subset(adiposity_proteins, adiposity == "WHR" & group == "Female")))
  
)

# protein-adiposity 
proteins_adiposity <- read.table("analysis/005_proteins_adiposity/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
proteins_adiposity <- proteins_adiposity[,c("sequence_id", "protein", "gene", "outcome", "group", "b")]
colnames(proteins_adiposity)[4] <- "adiposity"

table2 <- data.frame(
  analysis = c("adiposity-proteins > adiposity"),
  
  bmi_combined = c(nrow(subset(proteins_adiposity, adiposity == "BMI" & group == "Sex-combined"))),
  bmi_male = c(nrow(subset(proteins_adiposity, adiposity == "BMI" & group == "Male"))),
  bmi_female = c(nrow(subset(proteins_adiposity, adiposity == "BMI" & group == "Female"))),
  
  whr_combined = c(nrow(subset(proteins_adiposity, adiposity == "WHR" & group == "Sex-combined"))),
  whr_male = c(nrow(subset(proteins_adiposity, adiposity == "WHR" & group == "Male"))),
  whr_female = c(nrow(subset(proteins_adiposity, adiposity == "WHR" & group == "Female")))
)

## join
data <- left_join(proteins_adiposity, adiposity_proteins, by = c("sequence_id", "protein", "gene", "adiposity", "group"))
data <- data[complete.cases(data), ]
exclude_bmi <- subset(data, adiposity == "BMI")
exclude_bmi_combined <- subset(exclude_bmi, group == "Sex-combined")
exclude_bmi_combined <- exclude_bmi_combined$sequence_id
exclude_bmi_male <- subset(exclude_bmi, group == "Male")
exclude_bmi_male <- exclude_bmi_male$sequence_id
exclude_bmi_female <- subset(exclude_bmi, group == "Female")
exclude_bmi_female <- exclude_bmi_female$sequence_id

exclude_whr <- subset(data, adiposity == "WHR")
exclude_whr_combined <- subset(exclude_whr, group == "Sex-combined")
exclude_whr_combined <- exclude_whr_combined$sequence_id
exclude_whr_male <- subset(exclude_whr, group == "Male")
exclude_whr_male <- exclude_whr_male$sequence_id
exclude_whr_female <- subset(exclude_whr, group == "Female")
exclude_whr_female <- exclude_whr_female$sequence_id

## make association lists 
bmi <- subset(adiposity_proteins, adiposity == "BMI")
bmi_combined <- subset(bmi, group == "Sex-combined")
bmi_male <- subset(bmi, group == "Male")
bmi_female <- subset(bmi, group == "Female")

whr <- subset(adiposity_proteins, adiposity == "WHR")
whr_combined <- subset(whr, group == "Sex-combined")
whr_male <- subset(whr, group == "Male")
whr_female <- subset(whr, group == "Female")

## remove reverse associations
bmi_combined <- bmi_combined[!bmi_combined$sequence_id %in% exclude_bmi_combined , ]
bmi_male <- bmi_male[!bmi_male$sequence_id %in% exclude_bmi_male , ]
bmi_female <- bmi_female[!bmi_female$sequence_id %in% exclude_bmi_female , ]
whr_combined <- whr_combined[!whr_combined$sequence_id %in% exclude_whr_combined , ]
whr_male <- whr_male[!whr_male$sequence_id %in% exclude_whr_male , ]
whr_female <- whr_female[!whr_female$sequence_id %in% exclude_whr_female , ]

table3 <- data.frame(
  analysis = c("adiposity-proteins conflict"),
  
  bmi_combined = c(table1[1,2] - nrow(bmi_combined)),
  bmi_male = c(table1[1,3] - nrow(bmi_male)),
  bmi_female = c(table1[1,4] - nrow(bmi_female)),
  
  whr_combined = c(table1[1,5] - nrow(whr_combined)),
  whr_male = c(table1[1,6] - nrow(whr_male)),
  whr_female = c(table1[1,7] - nrow(whr_female))
)

table4 <- data.frame(
  analysis = c("adiposity-protein association"),
  
  bmi_combined = c(nrow(bmi_combined)),
  bmi_male = c(nrow(bmi_male)),
  bmi_female = c(nrow(bmi_female)),
  
  whr_combined = c(nrow(whr_combined)),
  whr_male = c(nrow(whr_male)),
  whr_female = c(nrow(whr_female))
)

## combine and save
adiposity_protein <- bind_rows(bmi_combined, bmi_female, bmi_male,
                  whr_combined, whr_female, whr_male)
write.table(adiposity_protein, "analysis/007_adiposity_proteins_colorectal/adiposity_protein.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# protein-cancer associations ====
## protein-cancer
protein_cancer <- read.table("analysis/011_protein_colorectal_cissnp/001_MR_results.txt", header = T, sep = "\t")
protein_cancer <- subset(protein_cancer, outcome != "overall-HRC-EUR")
protein_cancer <- subset(protein_cancer, outcome != "early_onset")
length(unique(protein_cancer$sequence_id))
table5 <- data.frame(table(protein_cancer$group, protein_cancer$outcome))
protein_cancer_id <- paste0(protein_cancer$sequence_id, "_", protein_cancer$outcome, "_", protein_cancer$group)
protein_cancer_id <- unique(protein_cancer_id)
protein_cancer <- subset(protein_cancer, pval < 0.05/1293)
protein_cancer <- protein_cancer[,c("sequence_id", "protein", "gene", "group", "outcome", "b")]
protein_cancer <- subset(protein_cancer, outcome != "overall-HRC-EUR")
protein_cancer <- subset(protein_cancer, outcome != "early_onset")
protein_cancer$id <- paste0(protein_cancer$sequence_id, "_", protein_cancer$outcome, "_", protein_cancer$group)
protein_cancer_associations <- paste0(protein_cancer$sequence_id, "_", protein_cancer$outcome, "_", protein_cancer$group)
table5 <- data.frame(table(protein_cancer$outcome, protein_cancer$group))

## cancer-protein 
cancer_protein <- read.table("analysis/006_colorectal_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
cancer_protein <- cancer_protein[,c("sequence_id", "protein", "gene", "group", "exposure", "b")]
cancer_protein <- subset(cancer_protein, exposure != "overall-HRC-EUR")
cancer_protein <- subset(cancer_protein, exposure != "early_onset")
cancer_protein$id <- paste0(cancer_protein$sequence_id, "_", cancer_protein$exposure, "_", cancer_protein$group)
cancer_protein_associations <- paste0(cancer_protein$sequence_id, "_", cancer_protein$exposure, "_", cancer_protein$group)
cancer_protein <- cancer_protein[cancer_protein$id %in% protein_cancer_id, ]
table6 <- data.frame(table(cancer_protein$exposure, cancer_protein$group))

## combine with coloc and save
protein_cancer <- protein_cancer[!protein_cancer$id %in% cancer_protein_associations, ]
protein_cancer$id <- paste0(protein_cancer$sequence_id, "_", protein_cancer$gene, "_", protein_cancer$protein, "_", protein_cancer$outcome, "_", protein_cancer$group)

colocalization <- read.table("analysis/009_colocalisation/results/002_coloc_results.txt", header = T, sep = "\t")
colocalization <- subset(colocalization , h4 > 0.8)
colocalization_associations <- colocalization$id
a <- protein_cancer[!protein_cancer$id %in% colocalization_associations, ]
table(a$outcome, a$group)

protein_cancer <- protein_cancer[protein_cancer$id %in% colocalization_associations, ]
write.table(protein_cancer, "analysis/007_adiposity_proteins_colorectal/protein_cancer.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# triplets ====
adiposity_protein <- select(adiposity_protein, adiposity, sequence_id, protein, gene, group, b)
protein_cancer <- select(protein_cancer, outcome, sequence_id, group, b)
data <- inner_join(protein_cancer, adiposity_protein, by = c("sequence_id", "group"))
table(data$group, data$outcome)
data_intermediate <- data[data$b.x * data$b.y >= 0, ]
write.table(data_intermediate, "analysis/007_adiposity_proteins_colorectal/adiposity_protein_cancer_consistent.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- select(data, adiposity, sequence_id, protein, gene, outcome, group, b.y, b.x)
colnames(data) <- c("adiposity", "sequence_id", "protein", "gene", "outcome", "group", "b_adiposity-protein", "b_protein-cancer")
data$id <- paste0(data$sequence_id, "_", data$gene, "_", data$protein)
write.table(data, "analysis/007_adiposity_proteins_colorectal/adiposity_protein_cancer_all.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
