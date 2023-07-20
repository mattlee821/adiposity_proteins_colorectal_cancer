rm(list=ls())
## set environment ====
library(dplyr)
library(tidyr)
library(functions)

# data ====
data <- read.table("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T, sep = "\t")

# make exposure_method id column
data$id <- paste0(data$exposure, "_", data$group, "_", data$method)

# make dataframe
data <- data[, c("id", "outcome", "b", "group", "exposure")]

# split into sex and exposure
data_male <- subset(data, group == "Male")
data_bmi_male <- subset(data_male, exposure == "BMI")
data_bmi_male <- data_bmi_male[,c("id", "outcome", "b")]
data_whr_male <- subset(data_male, exposure == "WHR")
data_whr_male <- data_whr_male[,c("id", "outcome", "b")]

data_female <- subset(data, group == "Female")
data_bmi_female <- subset(data_female, exposure == "BMI")
data_bmi_female <- data_bmi_female[,c("id", "outcome", "b")]
data_whr_female <- subset(data_female, exposure == "WHR")
data_whr_female <- data_whr_female[,c("id", "outcome", "b")]

data_combined <- subset(data, group == "Sex-combined")
data_bmi_combined <- subset(data_combined, exposure == "BMI")
data_bmi_combined <- data_bmi_combined[,c("id", "outcome", "b")]
data_whr_combined <- subset(data_combined, exposure == "WHR")
data_whr_combined <- data_whr_combined[,c("id", "outcome", "b")]

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
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Male"
data_bmi_male <- data1
data_bmi_male <- subset(data_bmi_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_bmi_female)
rownames(data1) <- data1[,1]
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Female"
data_bmi_female <- data1
data_bmi_female <- subset(data_bmi_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_bmi_combined)
rownames(data1) <- data1[,1]
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "BMI"
data1$group <- "Sex-combined"
data_bmi_combined <- data1
data_bmi_combined <- subset(data_bmi_combined, direction_group != "inconsistent")

## whr ====
### male
data1 <- as.data.frame(data_whr_male)
rownames(data1) <- data1[,1]
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Male"
data_whr_male <- data1
data_whr_male <- subset(data_whr_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_whr_female)
rownames(data1) <- data1[,1]
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Female"
data_whr_female <- data1
data_whr_female <- subset(data_whr_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_whr_combined)
rownames(data1) <- data1[,1]
data1 <- directional_consistency(data1)
data1 <- tibble::rownames_to_column(data1, "outcome")
data1 <- data1[,c("outcome", "direction_group")]
data1$exposure <- "whr"
data1$group <- "Sex-combined"
data_whr_combined <- data1
data_whr_combined <- subset(data_whr_combined, direction_group != "inconsistent")

# save directionally associated protein lists ====
a <- bind_rows(data_bmi_combined, data_bmi_female, data_bmi_male,
               data_whr_combined, data_whr_female, data_whr_male)
write.table(a, "analysis/001_adiposity_colorectal/001_directional_consistency.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
