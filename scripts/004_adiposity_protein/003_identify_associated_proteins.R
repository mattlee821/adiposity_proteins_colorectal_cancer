rm(list=ls())
## set environment ====
library(dplyr)
library(tidyr)

# data ====
data <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")

# make exposure_method id column
data$id <- paste0(data$exposure, "_", data$group, "_", data$method)

# make dataframe
data <- data[, c("id", "sequence_id", "b", "group", "exposure")]

# split into sex and exposure
data_male <- subset(data, group == "Male")
data_bmi_male <- subset(data_male, exposure == "BMI")
data_bmi_male <- data_bmi_male[,c("id", "sequence_id", "b")]
data_whr_male <- subset(data_male, exposure == "WHR")
data_whr_male <- data_whr_male[,c("id", "sequence_id", "b")]

data_female <- subset(data, group == "Female")
data_bmi_female <- subset(data_female, exposure == "BMI")
data_bmi_female <- data_bmi_female[,c("id", "sequence_id", "b")]
data_whr_female <- subset(data_female, exposure == "WHR")
data_whr_female <- data_whr_female[,c("id", "sequence_id", "b")]

data_combined <- subset(data, group == "Sex-combined")
data_bmi_combined <- subset(data_combined, exposure == "BMI")
data_bmi_combined <- data_bmi_combined[,c("id", "sequence_id", "b")]
data_whr_combined <- subset(data_combined, exposure == "WHR")
data_whr_combined <- data_whr_combined[,c("id", "sequence_id", "b")]

# convert to wide
data_bmi_male <- pivot_wider(data_bmi_male, names_from = id, values_from = b)
data_whr_male <- pivot_wider(data_whr_male, names_from = id, values_from = b)

data_bmi_female <- pivot_wider(data_bmi_female, names_from = id, values_from = b)
data_whr_female <- pivot_wider(data_whr_female, names_from = id, values_from = b)

data_bmi_combined <- pivot_wider(data_bmi_combined, names_from = id, values_from = b)
data_whr_combined <- pivot_wider(data_whr_combined, names_from = id, values_from = b)

# assign directions ====
## bmi ====
### male
data1 <- as.data.frame(data_bmi_male)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Male"
data_bmi_male <- data1
data_bmi_male <- subset(data_bmi_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_bmi_female)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Female"
data_bmi_female <- data1
data_bmi_female <- subset(data_bmi_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_bmi_combined)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Sex-combined"
data_bmi_combined <- data1
data_bmi_combined <- subset(data_bmi_combined, direction_group != "inconsistent")

## whr ====
### male
data1 <- as.data.frame(data_whr_male)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Male"
data_whr_male <- data1
data_whr_male <- subset(data_whr_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_whr_female)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Female"
data_whr_female <- data1
data_whr_female <- subset(data_whr_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_whr_combined)
rownames(data1) <- data1[,1]
data1 <- data1[,c(2:ncol(data1))]
data1$direction <- sapply(1:nrow(data1), function(x) ifelse(all(sign(data1[x,]) == 1), 1, ifelse(all(sign(data1[x,]) == -1), 2, 0)))
data1 <- data1 %>%
  mutate(direction_group = case_when(direction == 1 ~ "positive",
                               direction == 2 ~ "negative",
                               direction == 0 ~ "inconsistent"))
data1$direction_group <- factor(data1$direction_group,
                            levels = c("positive", "negative", "inconsistent"),ordered = TRUE)
data1 <- tibble::rownames_to_column(data1, "sequence_id")
data1 <- data1[,c("sequence_id", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Sex-combined"
data_whr_combined <- data1
data_whr_combined <- subset(data_whr_combined, direction_group != "inconsistent")

# save directionally associated protein lists ====
a <- subset(data_bmi_male, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_bmi_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_whr_male, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_whr_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_bmi_female, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_bmi_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_whr_female, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_whr_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_bmi_combined, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_bmi_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_whr_combined, direction_group != "inconsistent")
write.table(a, "analysis/002_adiposity_proteins/002_whr_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# pvalue thresholds ====
# data ====
data <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")

# make exposure_method id column
data$id <- paste0(data$exposure, "_", data$group, "_", data$method)

# make dataframe ====
data <- subset(data, method == "IVW-MRE")
data <- data[, c("id", "sequence_id", "pval", "group", "exposure")]

## subset for p 0.05 and phenospd
data1 <- subset(data, pval < 0.05)
data2 <- subset(data, pval < 0.05/1293)

# split into sex and exposure ====
## p < 0.05
data_male <- subset(data1, group == "Male")
data_bmi_male <- subset(data_male, exposure == "BMI")
data_whr_male <- subset(data_male, exposure == "WHR")

data_female <- subset(data1, group == "Female")
data_bmi_female <- subset(data_female, exposure == "BMI")
data_whr_female <- subset(data_female, exposure == "WHR")

data_combined <- subset(data1, group == "Sex-combined")
data_bmi_combined <- subset(data_combined, exposure == "BMI")
data_whr_combined <- subset(data_combined, exposure == "WHR")

table_p1 <- data.frame(
    exposure = c("bmi", "whr", "unique"),
    male_nominal = c(nrow(data_bmi_male), nrow(data_whr_male), length(unique(data_male$sequence_id))),
    female_nominal = c(nrow(data_bmi_female), nrow(data_whr_female), length(unique(data_female$sequence_id))),
    combined_nominal = c(nrow(data_bmi_combined), nrow(data_whr_combined), length(unique(data_combined$sequence_id)))
    )

## p < 0.05/1293
data_male <- subset(data2, group == "Male")
data_bmi_male <- subset(data_male, exposure == "BMI")
data_whr_male <- subset(data_male, exposure == "WHR")

data_female <- subset(data2, group == "Female")
data_bmi_female <- subset(data_female, exposure == "BMI")
data_whr_female <- subset(data_female, exposure == "WHR")

data_combined <- subset(data2, group == "Sex-combined")
data_bmi_combined <- subset(data_combined, exposure == "BMI")
data_whr_combined <- subset(data_combined, exposure == "WHR")

table_p2 <- data.frame(
    exposure = c("bmi", "whr", "unique"),
    male_phenospd = c(nrow(data_bmi_male), nrow(data_whr_male), length(unique(data_male$sequence_id))),
    female_phenospd = c(nrow(data_bmi_female), nrow(data_whr_female), length(unique(data_female$sequence_id))),
    combined_phenospd = c(nrow(data_bmi_combined), nrow(data_whr_combined), length(unique(data_combined$sequence_id)))
    )

## combine p tables 
table <- left_join(table_p1, table_p2)
write.table(table, "analysis/002_adiposity_proteins/table_N_p_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# associated directionally and p ====
data <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")
data <- subset(data, method == "IVW-MRE")
data <- subset(data, pval < 0.05/1293)

# split into sex and exposure
data_male <- subset(data, group == "Male")
data_bmi_male <- subset(data_male, exposure == "BMI")
data_whr_male <- subset(data_male, exposure == "WHR")

data_female <- subset(data, group == "Female")
data_bmi_female <- subset(data_female, exposure == "BMI")
data_whr_female <- subset(data_female, exposure == "WHR")

data_combined <- subset(data, group == "Sex-combined")
data_bmi_combined <- subset(data_combined, exposure == "BMI")
data_whr_combined <- subset(data_combined, exposure == "WHR")

# load directionally associated dataframes 
directionally_associated_bmi_male <- read.table("analysis/002_adiposity_proteins/002_bmi_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_bmi_male <- directionally_associated_bmi_male$sequence_id
directionally_associated_whr_male <- read.table("analysis/002_adiposity_proteins/002_whr_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_whr_male <- directionally_associated_whr_male$sequence_id
directionally_associated_bmi_female <- read.table("analysis/002_adiposity_proteins/002_bmi_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_bmi_female <- directionally_associated_bmi_female$sequence_id
directionally_associated_whr_female <- read.table("analysis/002_adiposity_proteins/002_whr_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_whr_female <- directionally_associated_whr_female$sequence_id
directionally_associated_bmi_combined <- read.table("analysis/002_adiposity_proteins/002_bmi_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_bmi_combined <- directionally_associated_bmi_combined$sequence_id
directionally_associated_whr_combined <- read.table("analysis/002_adiposity_proteins/002_whr_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_whr_combined <- directionally_associated_whr_combined$sequence_id

# pull out directionally consistent from p value a data frames
data_bmi_male <- data_bmi_male[data_bmi_male$sequence_id %in% directionally_associated_bmi_male, ]
data_whr_male <- data_whr_male[data_whr_male$sequence_id %in% directionally_associated_whr_male, ]

data_bmi_female <- data_bmi_female[data_bmi_female$sequence_id %in% directionally_associated_bmi_female, ]
data_whr_female <- data_whr_female[data_whr_female$sequence_id %in% directionally_associated_whr_female, ]

data_bmi_combined <- data_bmi_combined[data_bmi_combined$sequence_id %in% directionally_associated_bmi_combined, ]
data_whr_combined <- data_whr_combined[data_whr_combined$sequence_id %in% directionally_associated_whr_combined, ]

# make table of N associations
table_p3 <- data.frame(
    exposure = c("bmi", "whr", "unique"),
    male_phenospd_consistent = c(nrow(data_bmi_male), nrow(data_whr_male), length(unique(data_male$sequence_id))),
    female_phenospd_consistent = c(nrow(data_bmi_female), nrow(data_whr_female), length(unique(data_female$sequence_id))),
    combined_phenospd_consistent = c(nrow(data_bmi_combined), nrow(data_whr_combined), length(unique(data_combined$sequence_id)))
    )

table <- left_join(table, table_p3)
write.table(table, "analysis/002_adiposity_proteins/table_N_p_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make table of associations ====
data_associations <- bind_rows(
    data_bmi_male, data_whr_male, 
    data_bmi_female, data_whr_female, 
    data_bmi_combined, data_whr_combined)

write.table(data_associations, "analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# unqiue associations table
data <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")

male <- subset(data, group == "Male")
female <- subset(data, group == "Female")
combined <- subset(data, group == "Sex-combined")

male_proteins <- unique(male$protein)
female_proteins <- unique(female$protein)
combined_proteins <- unique(combined$protein)

male_female_proteins <- c(male_proteins, female_proteins)
male_combined_proteins <- c(male_proteins, combined_proteins)
female_combined_proteins <- c(female_proteins, combined_proteins)

male_not_female <- male_proteins[!male_proteins %in% female_proteins]
male_not_combined <- male_proteins[!male_proteins %in% combined_proteins]
male_unique <- male_proteins[!male_proteins %in% female_combined_proteins]

female_not_male <- female_proteins[!female_proteins %in% male_proteins]
female_not_combined <- female_proteins[!female_proteins %in% combined_proteins]
female_unique <- female_proteins[!female_proteins %in% male_combined_proteins]

combined_not_male <- combined_proteins[!combined_proteins %in% male_proteins]
combined_not_female <- combined_proteins[!combined_proteins %in% female_proteins]
combined_unique <- combined_proteins[!combined_proteins %in% male_female_proteins]

unique_proteins <- c(male_proteins, female_proteins, combined_proteins)
unique_proteins <- unique(unique_proteins)

table <- data.frame(
  row.names = c("N", "male", "female", "combined", "unique"),
  N = c(4907, length(male_proteins), length(female_proteins), length(combined_proteins),length(unique_proteins)),
  male = c(length(male_proteins),NA,length(male_not_female),length(male_not_combined),length(male_unique)),
  female = c(length(female_proteins),length(female_not_male),NA,length(female_not_combined),length(female_unique)),
  combined = c(length(combined_proteins),length(combined_not_male),length(combined_not_female),NA,length(combined_unique))
)

write.table(table, "analysis/002_adiposity_proteins/table_N_unique_associations.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


