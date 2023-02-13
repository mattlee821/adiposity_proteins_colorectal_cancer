rm(list=ls())
library(dplyr)

# adiposity protein associations ====
adiposity_proteins <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
adiposity_proteins <- adiposity_proteins[,c("sequence_id", "protein", "gene", "exposure", "group")]
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

# protein-adiposity associations ====
proteins_adiposity <- read.table("analysis/005_proteins_adiposity/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
proteins_adiposity <- proteins_adiposity[,c("sequence_id", "protein", "gene", "outcome", "group")]
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

# join ====
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

# make association lists ====
bmi <- subset(adiposity_proteins, adiposity == "BMI")
bmi_combined <- subset(bmi, group == "Sex-combined")
bmi_male <- subset(bmi, group == "Male")
bmi_female <- subset(bmi, group == "Female")

whr <- subset(adiposity_proteins, adiposity == "WHR")
whr_combined <- subset(whr, group == "Sex-combined")
whr_male <- subset(whr, group == "Male")
whr_female <- subset(whr, group == "Female")

# remove reverse associations
## BMI
bmi_combined <- bmi_combined[!bmi_combined$sequence_id %in% exclude_bmi_combined , ]
write.table(bmi_combined, "analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

bmi_male <- bmi_male[!bmi_male$sequence_id %in% exclude_bmi_male , ]
write.table(bmi_male, "analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

bmi_female <- bmi_female[!bmi_female$sequence_id %in% exclude_bmi_female , ]
write.table(bmi_female, "analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## whr
whr_combined <- whr_combined[!whr_combined$sequence_id %in% exclude_whr_combined , ]
write.table(whr_combined, "analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

whr_male <- whr_male[!whr_male$sequence_id %in% exclude_whr_male , ]
write.table(whr_male, "analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

whr_female <- whr_female[!whr_female$sequence_id %in% exclude_whr_female , ]
write.table(whr_female, "analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

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

# adiposity-protein cancer associations ====
protein_cancer <- read.table("analysis/003_proteins_colorectal/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
protein_cancer <- protein_cancer[,c("sequence_id", "protein", "gene", "group", "outcome")]

# make  association list bmi ====
data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_combined.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
a <- subset(a, outcome != "early_onset")
bmi_protein_cancer_combined <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_bmi_protein_crc_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_male.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
bmi_protein_cancer_male <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_bmi_protein_crc_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_bmi_proteins_female.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
bmi_protein_cancer_female <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_bmi_protein_crc_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make  association list whr ====
data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_combined.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
whr_protein_cancer_combined <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_whr_protein_crc_combined.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_male.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
whr_protein_cancer_male <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_whr_protein_crc_male.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/000_whr_proteins_female.txt", header = T, sep = "\t")
a <- left_join(data, protein_cancer, by = c("sequence_id", "protein", "gene", "group"))
a <- a[complete.cases(a), ]
a <- subset(a, outcome != "distal")
whr_protein_cancer_female <- unique(a$sequence_id)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/001_whr_protein_crc_female.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# table ====
table5 <- data.frame(
  analysis = c("adiposity-proteins > cancer"),
  
  bmi_combined = c(length(bmi_protein_cancer_combined)),
  bmi_male = c(length(bmi_protein_cancer_male)),
  bmi_female = c(length(bmi_protein_cancer_female)),
  
  whr_combined = c(length(whr_protein_cancer_combined)),
  whr_male = c(length(whr_protein_cancer_male)),
  whr_female = c(length(whr_protein_cancer_female))
  
)

# remove crc-associated proteins from list and make final association table ====
filenames <- dir("analysis/007_adiposity_proteins_colorectal/phenospd/", recursive=TRUE, full.names=TRUE, pattern="001")
list <- lapply(filenames, read.table, header = T, sep = "\t") 
adiposity_proteins_crc <- bind_rows(list)
colnames(adiposity_proteins_crc)[6] <- "cancer"

crc_proteins <- read.table("analysis/006_colorectal_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
colnames(crc_proteins)[1] <- "cancer"
crc_proteins <- crc_proteins[,c("sequence_id", "cancer", "group")]
crc_proteins$reverse <- "yes"

data <- left_join(adiposity_proteins_crc, crc_proteins, by = c("sequence_id", "cancer", "group"))

# table 
table6 <- data.frame(
  analysis = c("adiposity-protein < cancer"),
  
  bmi_combined = c(nrow(subset(data, adiposity == "BMI" & group == "Sex-combined" & reverse == "yes"))),
  bmi_male = c(nrow(subset(data, adiposity == "BMI" & group == "Male" & reverse == "yes"))),
  bmi_female = c(nrow(subset(data, adiposity == "BMI" & group == "Female" & reverse == "yes"))),
  
  whr_combined = c(nrow(subset(data, adiposity == "whr" & group == "Sex-combined" & reverse == "yes"))),
  whr_male = c(nrow(subset(data, adiposity == "whr" & group == "Male" & reverse == "yes"))),
  whr_female = c(nrow(subset(data, adiposity == "whr" & group == "Female" & reverse == "yes")))
)

# table 
data_forward <- data[is.na(data$reverse),]
nrow(subset(data_forward, adiposity == "BMI" & group == "Male"))

table7 <- data.frame(
  analysis = c("adiposity-protein-cancer associations"),
  
  bmi_combined = nrow(subset(data_forward, adiposity == "BMI" & group == "Sex-combined")),
  bmi_male = nrow(subset(data_forward, adiposity == "BMI" & group == "Male")),
  bmi_female = nrow(subset(data_forward, adiposity == "BMI" & group == "Female")),
  
  whr_combined = nrow(subset(data_forward, adiposity == "WHR" & group == "Sex-combined")),
  whr_male = nrow(subset(data_forward, adiposity == "WHR" & group == "Male")),
  whr_female = nrow(subset(data_forward, adiposity == "WHR" & group == "Female"))
)

table7a <- data.frame(
  analysis = c("adiposity-protein-cancer associations - N unique proteins"),
  
  bmi_combined = c(table5[1,2] - table6[1,2]),
  bmi_male = c(table5[1,3] - table6[1,3]),
  bmi_female = c(table5[1,4] - table6[1,4]),
  
  whr_combined = c(table5[1,5] - table6[1,5]),
  whr_male = c(table5[1,6] - table6[1,6]),
  whr_female = c(table5[1,7] - table6[1,7])
)

# save
data <- data[!complete.cases(data), ]
data <- data[,c("adiposity", "sequence_id", "protein", "gene", "cancer", "group")]
write.table(data, "analysis/007_adiposity_proteins_colorectal/phenospd/002_adiposity_protein_crc_associations.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
adiposity_protein_cancer <- data

b <- subset(data, group == "Sex-combined" & adiposity == "BMI")
nrow(b)
length(unique(b$sequence_id))
length(unique(b$cancer))

# make pos/neg list for MVMR ====
data <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/002_adiposity_protein_crc_associations.txt", header = T, sep = "\t")

## make a list of positive associations (increase-increase) ====
adiposity_protein <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
adiposity_protein <- adiposity_protein[,c("exposure", "sequence_id", "protein", "gene", "group", "b")]
colnames(adiposity_protein)[1] <- "adiposity"
colnames(adiposity_protein)[6] <- "b_adiposity_protein"

adiposity_protein_pos <- subset(adiposity_protein, b_adiposity_protein > 0)
adiposity_protein_pos <- left_join(adiposity_protein_pos, data, by = c("adiposity", "sequence_id", "protein", "gene", "group"))
adiposity_protein_pos <- adiposity_protein_pos[complete.cases(adiposity_protein_pos), ]

protein_cancer <- read.table("analysis/003_proteins_colorectal/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
protein_cancer$OR <- exp(protein_cancer$b)
protein_cancer <- protein_cancer[,c("outcome", "sequence_id", "protein", "gene", "group", "OR")]
colnames(protein_cancer)[1] <- "cancer"
colnames(protein_cancer)[6] <- "OR_protein_cancer"

protein_cancer_pos <- subset(protein_cancer, OR_protein_cancer > 1)
protein_cancer_pos <- left_join(protein_cancer_pos, data, by = c("cancer", "sequence_id", "protein", "gene", "group"))
protein_cancer_pos <- protein_cancer_pos[complete.cases(protein_cancer_pos), ]

adiposity_protein_cancer_pos <- left_join(adiposity_protein_pos, protein_cancer_pos, by = c("adiposity", "cancer", "sequence_id", "protein", "gene", "group"))
adiposity_protein_cancer_pos <- adiposity_protein_cancer_pos[complete.cases(adiposity_protein_cancer_pos), ]
a <- unique(adiposity_protein_cancer_pos$sequence_id)

write.table(adiposity_protein_cancer_pos, "analysis/007_adiposity_proteins_colorectal/phenospd/003_adiposity_protein_cancer_pos.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## make a list of negative associations (decrease-decrease) ====
adiposity_protein <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
adiposity_protein <- adiposity_protein[,c("exposure", "sequence_id", "protein", "gene", "group", "b")]
colnames(adiposity_protein)[1] <- "adiposity"
colnames(adiposity_protein)[6] <- "b_adiposity_protein"

adiposity_protein_neg <- subset(adiposity_protein, b_adiposity_protein < 0)
adiposity_protein_neg <- left_join(adiposity_protein_neg, data, by = c("adiposity", "sequence_id", "protein", "gene", "group"))
adiposity_protein_neg <- adiposity_protein_neg[complete.cases(adiposity_protein_neg), ]

protein_cancer <- read.table("analysis/003_proteins_colorectal/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
protein_cancer$OR <- exp(protein_cancer$b)
protein_cancer <- protein_cancer[,c("outcome", "sequence_id", "protein", "gene", "group", "OR")]
colnames(protein_cancer)[1] <- "cancer"
colnames(protein_cancer)[6] <- "OR_protein_cancer"

protein_cancer_neg <- subset(protein_cancer, OR_protein_cancer < 1)
protein_cancer_neg <- left_join(protein_cancer_neg, data, by = c("cancer", "sequence_id", "protein", "gene", "group"))
protein_cancer_neg <- protein_cancer_neg[complete.cases(protein_cancer_neg), ]

adiposity_protein_cancer_neg <- left_join(adiposity_protein_neg, protein_cancer_neg, by = c("adiposity", "cancer", "sequence_id", "protein", "gene", "group"))
adiposity_protein_cancer_neg <- adiposity_protein_cancer_neg[complete.cases(adiposity_protein_cancer_neg), ]
a <- unique(adiposity_protein_cancer_neg$sequence_id)

write.table(adiposity_protein_cancer_neg, "analysis/007_adiposity_proteins_colorectal/phenospd/003_adiposity_protein_cancer_neg.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

nrow(adiposity_protein_cancer)
length(unique(adiposity_protein_cancer$sequence_id))
length(unique(adiposity_protein_cancer$cancer))

## combine ====
a <- rbind(adiposity_protein_cancer_neg, adiposity_protein_cancer_pos)
write.table(a, "analysis/007_adiposity_proteins_colorectal/phenospd/004_adiposity_protein_cancer_pos_neg.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

table <- select(a, adiposity, sequence_id, protein, gene, cancer, group)
write.table(table, "analysis/tables/007_adiposity_proteins_colorectal_intermediates.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

nrow(a)
length(unique(a$sequence_id))
length(unique(a$cancer))

nrow(adiposity_protein_cancer_pos)
length(unique(adiposity_protein_cancer_pos$sequence_id))
length(unique(adiposity_protein_cancer_pos$cancer))

nrow(adiposity_protein_cancer_neg)
length(unique(adiposity_protein_cancer_neg$sequence_id))
length(unique(adiposity_protein_cancer_neg$cancer))

# table ====
table8 <- data.frame(
  analysis = c("adiposity-protein-cancer intermediates"),
  
  bmi_combined = c(nrow(subset(a, adiposity == "BMI" & group == "Sex-combined"))),
  bmi_male = c(nrow(subset(a, adiposity == "BMI" & group == "Male"))),
  bmi_female = c(nrow(subset(a, adiposity == "BMI" & group == "Female"))),
  
  whr_combined = c(nrow(subset(a, adiposity == "WHR" & group == "Sex-combined"))),
  whr_male = c(nrow(subset(a, adiposity == "WHR" & group == "Male"))),
  whr_female = c(nrow(subset(a, adiposity == "WHR" & group == "Female")))
  
)

# unique intermediate proteins
bmi_combined <- subset(a, adiposity == "BMI" & group == "Sex-combined")
bmi_male <- subset(a, adiposity == "BMI" & group == "Male")
bmi_female <- subset(a, adiposity == "BMI" & group == "Female")

whr_combined <- subset(a, adiposity == "WHR" & group == "Sex-combined")
whr_male <- subset(a, adiposity == "WHR" & group == "Male")
whr_female <- subset(a, adiposity == "WHR" & group == "Female")

table8a <- data.frame(
  analysis = c("adiposity-protein-cancer intermediates - N unique proteins"),
  
  bmi_combined = length(unique(bmi_combined$sequence_id)),
  bmi_male = length(unique(bmi_male$sequence_id)),
  bmi_female = length(unique(bmi_female$sequence_id)),
  
  whr_combined = length(unique(whr_combined$sequence_id)),
  whr_male = length(unique(whr_male$sequence_id)),
  whr_female = length(unique(whr_female$sequence_id))
  
)

# MVMR by cancer
table9 <- data.frame(
  analysis = c("overall",
               "overall-EUR",
               "colon",
               "distal",
               "proximal",
               "rectal"),
  
  bmi_combined = c(
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "overall")),
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "colon")),
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "distal")),
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "proximal")),
    nrow(subset(a, adiposity == "BMI" & group == "Sex-combined" & cancer == "rectal"))
  ),
  
  bmi_male = c(  
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "overall")),
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "colon")),
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "distal")),
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "proximal")),
    nrow(subset(a, adiposity == "BMI" & group == "Male" & cancer == "rectal"))
  ),
  
  bmi_female = c(
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "overall")),
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "colon")),
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "distal")),
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "proximal")),
    nrow(subset(a, adiposity == "BMI" & group == "Female" & cancer == "rectal"))
  ),
  
  whr_combined = c(
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "overall")),
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "colon")),
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "distal")),
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "proximal")),
    nrow(subset(a, adiposity == "WHR" & group == "Sex-combined" & cancer == "rectal"))
  ),
  
  whr_male = c(  
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "overall")),
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "colon")),
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "distal")),
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "proximal")),
    nrow(subset(a, adiposity == "WHR" & group == "Male" & cancer == "rectal"))
  ),
  
  whr_female = c(
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "overall")),
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "overall-HRC-EUR")),
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "colon")),
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "distal")),
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "proximal")),
    nrow(subset(a, adiposity == "WHR" & group == "Female" & cancer == "rectal"))
  )
)

table <- rbind(table1,table2,table3,table4,table5,table6,table7,table7a,table8,table8a,table9)
write.table(table, "analysis/007_adiposity_proteins_colorectal/table_N_associations.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

