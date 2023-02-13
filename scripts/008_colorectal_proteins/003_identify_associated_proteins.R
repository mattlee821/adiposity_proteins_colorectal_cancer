rm(list=ls())
## set environment ====
library(dplyr)
library(tidyr)
library(forcats)

# data ====
data <- read.table("analysis/006_colorectal_proteins/001_MR_results.txt", header = T, sep = "\t")

# make exposure_method id column
data$id <- paste0(data$exposure, "_", data$group, "_", data$method)

# make dataframe
data <- data[, c("id", "sequence_id", "b", "group", "exposure")]

# split into sex and exposure
data_male <- subset(data, group == "Male")
data_colon_male <- subset(data_male, exposure == "colon")
data_colon_male <- data_colon_male[,c("id", "sequence_id", "b")]
data_distal_male <- subset(data_male, exposure == "distal")
data_distal_male <- data_distal_male[,c("id", "sequence_id", "b")]
data_proximal_male <- subset(data_male, exposure == "proximal")
data_proximal_male <- data_proximal_male[,c("id", "sequence_id", "b")]
data_rectal_male <- subset(data_male, exposure == "rectal")
data_rectal_male <- data_colon_male[,c("id", "sequence_id", "b")]
data_overall_male <- subset(data_male, exposure == "overall")
data_overall_male <- data_overall_male[,c("id", "sequence_id", "b")]

data_female <- subset(data, group == "Female")
data_colon_female <- subset(data_female, exposure == "colon")
data_colon_female <- data_colon_female[,c("id", "sequence_id", "b")]
data_distal_female <- subset(data_female, exposure == "distal")
data_distal_female <- data_distal_female[,c("id", "sequence_id", "b")]
data_proximal_female <- subset(data_female, exposure == "proximal")
data_proximal_female <- data_proximal_female[,c("id", "sequence_id", "b")]
data_rectal_female <- subset(data_female, exposure == "rectal")
data_rectal_female <- data_colon_female[,c("id", "sequence_id", "b")]
data_overall_female <- subset(data_female, exposure == "overall")
data_overall_female <- data_overall_female[,c("id", "sequence_id", "b")]

data_combined <- subset(data, group == "Sex-combined")
data_colon_combined <- subset(data_combined, exposure == "colon")
data_colon_combined <- data_colon_combined[,c("id", "sequence_id", "b")]
data_distal_combined <- subset(data_combined, exposure == "distal")
data_distal_combined <- data_distal_combined[,c("id", "sequence_id", "b")]
data_proximal_combined <- subset(data_combined, exposure == "proximal")
data_proximal_combined <- data_proximal_combined[,c("id", "sequence_id", "b")]
data_rectal_combined <- subset(data_combined, exposure == "rectal")
data_rectal_combined <- data_colon_combined[,c("id", "sequence_id", "b")]
data_overall_combined <- subset(data_combined, exposure == "overall")
data_overall_combined <- data_overall_combined[,c("id", "sequence_id", "b")]
data_overall_HRC_combined <- subset(data_combined, exposure == "overall-HRC-EUR")
data_overall_HRC_combined <- data_overall_HRC_combined[,c("id", "sequence_id", "b")]
data_overall_early_onset <- subset(data_combined, exposure == "early_onset")
data_overall_early_onset <- data_overall_early_onset[,c("id", "sequence_id", "b")]

# convert to wide
data_colon_male <- pivot_wider(data_colon_male, names_from = id, values_from = b)
data_distal_male <- pivot_wider(data_distal_male, names_from = id, values_from = b)
data_proximal_male <- pivot_wider(data_proximal_male, names_from = id, values_from = b)
data_rectal_male <- pivot_wider(data_rectal_male, names_from = id, values_from = b)
data_overall_male <- pivot_wider(data_overall_male, names_from = id, values_from = b)

data_colon_female <- pivot_wider(data_colon_female, names_from = id, values_from = b)
data_distal_female <- pivot_wider(data_distal_female, names_from = id, values_from = b)
data_proximal_female <- pivot_wider(data_proximal_female, names_from = id, values_from = b)
data_rectal_female <- pivot_wider(data_rectal_female, names_from = id, values_from = b)
data_overall_female <- pivot_wider(data_overall_female, names_from = id, values_from = b)

data_colon_combined <- pivot_wider(data_colon_combined, names_from = id, values_from = b)
data_distal_combined <- pivot_wider(data_distal_combined, names_from = id, values_from = b)
data_proximal_combined <- pivot_wider(data_proximal_combined, names_from = id, values_from = b)
data_rectal_combined <- pivot_wider(data_rectal_combined, names_from = id, values_from = b)
data_overall_combined <- pivot_wider(data_overall_combined, names_from = id, values_from = b)
data_overall_HRC_combined <- pivot_wider(data_overall_HRC_combined, names_from = id, values_from = b)
data_overall_early_onset <- pivot_wider(data_overall_early_onset, names_from = id, values_from = b)

# assign directions ====
## colon ====
### male
data1 <- as.data.frame(data_colon_male)
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
data1$exposure <- "colon"
data1$group <- "Male"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_colon_male <- data1
data_colon_male <- subset(data_colon_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_colon_female)
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
data1$exposure <- "colon"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_colon_female <- data1
data_colon_female <- subset(data_colon_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_colon_combined)
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
data1$exposure <- "colon"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_colon_combined <- data1
data_colon_combined <- subset(data_colon_combined, direction_group != "inconsistent")

## distal ====
### male
data1 <- as.data.frame(data_distal_male)
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
data1$exposure <- "distal"
data1$group <- "Male"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_distal_male <- data1
data_distal_male <- subset(data_distal_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_distal_female)
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
data1$exposure <- "distal"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_distal_female <- data1
data_distal_female <- subset(data_distal_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_distal_combined)
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
data1$exposure <- "distal"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_distal_combined <- data1
data_distal_combined <- subset(data_distal_combined, direction_group != "inconsistent")

## rectal ====
### male
data1 <- as.data.frame(data_rectal_male)
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
data1$exposure <- "rectal"
data1$group <- "Male"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_rectal_male <- data1
data_rectal_male <- subset(data_rectal_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_rectal_female)
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
data1$exposure <- "rectal"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_rectal_female <- data1
data_rectal_female <- subset(data_rectal_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_rectal_combined)
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
data1$exposure <- "rectal"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_rectal_combined <- data1
data_rectal_combined <- subset(data_rectal_combined, direction_group != "inconsistent")

## proximal ====
### male
data1 <- as.data.frame(data_proximal_male)
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
data1$exposure <- "proximal"
data1$group <- "Male"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_proximal_male <- data1
data_proximal_male <- subset(data_proximal_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_proximal_female)
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
data1$exposure <- "proximal"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_proximal_female <- data1
data_proximal_female <- subset(data_proximal_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_proximal_combined)
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
data1$exposure <- "proximal"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_proximal_combined <- data1
data_proximal_combined <- subset(data_proximal_combined, direction_group != "inconsistent")

## overall ====
### male
data1 <- as.data.frame(data_overall_male)
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
data1$exposure <- "overall"
data1$group <- "Male"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_overall_male <- data1
data_overall_male <- subset(data_overall_male, direction_group != "inconsistent")

### female
data1 <- as.data.frame(data_overall_female)
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
data1$exposure <- "overall"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_overall_female <- data1
data_overall_female <- subset(data_overall_female, direction_group != "inconsistent")

### combined
data1 <- as.data.frame(data_overall_combined)
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
data1$exposure <- "overall"
data1$group <- "Female"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_overall_combined <- data1
data_overall_combined <- subset(data_overall_combined, direction_group != "inconsistent")

data1 <- as.data.frame(data_overall_HRC_combined)
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
data1$exposure <- "overall"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_overall_HRC_combined <- data1
data_overall_HRC_combined <- subset(data_overall_HRC_combined, direction_group != "inconsistent")

data1 <- as.data.frame(data_overall_early_onset)
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
data1$exposure <- "overall"
data1$group <- "Sex-combined"
data1$direction_group <- fct_explicit_na(data1$direction_group, "NA")
data_overall_early_onset <- data1
data_overall_early_onset <- subset(data_overall_early_onset, direction_group != "inconsistent")

# make table of directionally associated N ====
table <- data.frame(
  exposure = c("overall", "overall_HRC", "early_onset", "colon", "distal", "proximal", "rectal"),
  
  combined_positive = c(table(data_overall_combined$direction_group)[[1]],
                        table(data_overall_HRC_combined$direction_group)[[1]],
                        table(data_overall_early_onset$direction_group)[[1]],
                        table(data_colon_combined$direction_group)[[1]],
                        table(data_distal_combined$direction_group)[[1]],
                        table(data_proximal_combined$direction_group)[[1]],
                        table(data_rectal_combined$direction_group)[[1]]),
  combined_negative = c(table(data_overall_combined$direction_group)[[2]],
                        table(data_overall_HRC_combined$direction_group)[[2]],
                        table(data_overall_early_onset$direction_group)[[2]],
                        table(data_colon_combined$direction_group)[[2]],
                        table(data_distal_combined$direction_group)[[2]],
                        table(data_proximal_combined$direction_group)[[2]],
                        table(data_rectal_combined$direction_group)[[2]]),
  combined_one_method = c(table(data_overall_combined$direction_group)[[4]],
                          NA,
                          table(data_overall_early_onset$direction_group)[[4]],
                          table(data_colon_combined$direction_group)[[4]],
                          NA,
                          table(data_proximal_combined$direction_group)[[4]],
                          table(data_rectal_combined$direction_group)[[4]]),
  
  male_positive = c(table(data_overall_male$direction_group)[[1]],
                    NA,
                    NA,
                    table(data_colon_male$direction_group)[[1]],
                    table(data_distal_male$direction_group)[[1]],
                    table(data_proximal_male$direction_group)[[1]],
                    table(data_rectal_male$direction_group)[[1]]),
  male_negative = c(table(data_overall_male$direction_group)[[2]],
                    NA,
                    NA,
                    table(data_colon_male$direction_group)[[2]],
                    table(data_distal_male$direction_group)[[2]],
                    table(data_proximal_male$direction_group)[[2]],
                    table(data_rectal_male$direction_group)[[2]]),
  male_one_method = c(NA,
                      NA,
                      NA,
                      NA,
                      NA,
                      table(data_proximal_male$direction_group)[[4]],
                      NA),
  
  female_positive = c(table(data_overall_female$direction_group)[[1]],
                      NA,
                      NA,
                      table(data_colon_female$direction_group)[[1]],
                      table(data_distal_female$direction_group)[[1]],
                      table(data_proximal_female$direction_group)[[1]],
                      table(data_rectal_female$direction_group)[[1]]),
  female_negative = c(table(data_overall_female$direction_group)[[2]],
                      NA,
                      NA,
                      table(data_colon_female$direction_group)[[2]],
                      table(data_distal_female$direction_group)[[2]],
                      table(data_proximal_female$direction_group)[[2]],
                      table(data_rectal_female$direction_group)[[2]]),
  female_one_method = c(NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        table(data_proximal_female$direction_group)[[4]],
                        NA)
)

# 
write.table(table, "analysis/006_colorectal_proteins/table_N_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# save directionally associated protein lists ====
a <- subset(data_colon_male, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_colon_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_distal_male, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_distal_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_proximal_male, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_proximal_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_rectal_male, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_rectal_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_overall_male, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_overall_male_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

a <- subset(data_colon_female, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_colon_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_distal_female, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_distal_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_proximal_female, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_proximal_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_rectal_female, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_rectal_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_overall_female, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_overall_female_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

a <- subset(data_colon_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_colon_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_distal_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_distal_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_proximal_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_proximal_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_rectal_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_rectal_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_overall_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_overall_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_overall_HRC_combined, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_overall_HRC_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
a <- subset(data_overall_early_onset, direction_group != "inconsistent")
write.table(a, "analysis/006_colorectal_proteins/002_overall_early_onset_combined_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# pvalue thresholds ====
# data ====
data <- read.table("analysis/006_colorectal_proteins/001_MR_results.txt", header = T, sep = "\t")

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
# split into sex and exposure
data_male <- subset(data1, group == "Male")
data_colon_male <- subset(data_male, exposure == "colon")
data_distal_male <- subset(data_male, exposure == "distal")
data_proximal_male <- subset(data_male, exposure == "proximal")
data_rectal_male <- subset(data_male, exposure == "rectal")
data_overall_male <- subset(data_male, exposure == "overall")

data_female <- subset(data1, group == "Female")
data_colon_female <- subset(data_female, exposure == "colon")
data_distal_female <- subset(data_female, exposure == "distal")
data_proximal_female <- subset(data_female, exposure == "proximal")
data_rectal_female <- subset(data_female, exposure == "rectal")
data_overall_female <- subset(data_female, exposure == "overall")

data_combined <- subset(data1, group == "Sex-combined")
data_colon_combined <- subset(data_combined, exposure == "colon")
data_distal_combined <- subset(data_combined, exposure == "distal")
data_proximal_combined <- subset(data_combined, exposure == "proximal")
data_rectal_combined <- subset(data_combined, exposure == "rectal")
data_overall_combined <- subset(data_combined, exposure == "overall")
data_overall_HRC_combined <- subset(data_combined, exposure == "overall-HRC-EUR")
data_overall_early_onset <- subset(data_combined, exposure == "early_onset")

# make table of N 
table_p1 <- data.frame(
  exposure = c("overall", "overall_HRC", "early_onset", "colon", "distal", "proximal", "rectal", "unique"),
  male_nominal = c(nrow(data_overall_male), 
                   NA,
                   NA,
                   nrow(data_colon_male), 
                   nrow(data_distal_male), 
                   nrow(data_proximal_male), 
                   nrow(data_rectal_male), 
                   length(unique(data_male$sequence_id))),
  female_nominal = c(nrow(data_overall_female), 
                     NA,
                     NA,
                     nrow(data_colon_female), 
                     nrow(data_distal_female), 
                     nrow(data_proximal_female), 
                     nrow(data_rectal_female), 
                     length(unique(data_female$sequence_id))),
  combined_nominal = c(nrow(data_overall_combined), 
                       nrow(data_overall_HRC_combined),
                       nrow(data_overall_early_onset),
                       nrow(data_colon_combined), 
                       nrow(data_distal_combined), 
                       nrow(data_proximal_combined), 
                       nrow(data_rectal_combined), 
                       length(unique(data_combined$sequence_id)))
)

## p < 0.05/1293
data_male <- subset(data2, group == "Male")
data_colon_male <- subset(data_male, exposure == "colon")
data_distal_male <- subset(data_male, exposure == "distal")
data_proximal_male <- subset(data_male, exposure == "proximal")
data_rectal_male <- subset(data_male, exposure == "rectal")
data_overall_male <- subset(data_male, exposure == "overall")

data_female <- subset(data2, group == "Female")
data_colon_female <- subset(data_female, exposure == "colon")
data_distal_female <- subset(data_female, exposure == "distal")
data_proximal_female <- subset(data_female, exposure == "proximal")
data_rectal_female <- subset(data_female, exposure == "rectal")
data_overall_female <- subset(data_female, exposure == "overall")

data_combined <- subset(data2, group == "Sex-combined")
data_colon_combined <- subset(data_combined, exposure == "colon")
data_distal_combined <- subset(data_combined, exposure == "distal")
data_proximal_combined <- subset(data_combined, exposure == "proximal")
data_rectal_combined <- subset(data_combined, exposure == "rectal")
data_overall_combined <- subset(data_combined, exposure == "overall")
data_overall_HRC_combined <- subset(data_combined, exposure == "overall-HRC-EUR")
data_overall_early_onset <- subset(data_combined, exposure == "early_onset")

# make table of N 
table_p2 <- data.frame(
  exposure = c("overall", "overall_HRC", "early_onset", "colon", "distal", "proximal", "rectal", "unique"),
  male_phenospd = c(nrow(data_overall_male), 
                      NA,
                    NA,
                      nrow(data_colon_male), 
                      nrow(data_distal_male), 
                      nrow(data_proximal_male), 
                      nrow(data_rectal_male), 
                      length(unique(data_male$sequence_id))),
  female_phenospd = c(nrow(data_overall_female), 
                        NA,
                      NA,
                        nrow(data_colon_female), 
                        nrow(data_distal_female), 
                        nrow(data_proximal_female), 
                        nrow(data_rectal_female), 
                        length(unique(data_female$sequence_id))),
  combined_phenospd = c(nrow(data_overall_combined), 
                          nrow(data_overall_HRC_combined),
                        nrow(data_overall_early_onset),
                        nrow(data_colon_combined), 
                          nrow(data_distal_combined), 
                          nrow(data_proximal_combined), 
                          nrow(data_rectal_combined), 
                          length(unique(data_combined$sequence_id)))
)

## combine p tables 
table <- left_join(table_p1, table_p2)
write.table(table, "analysis/006_colorectal_proteins/table_N_p_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# associated directionally and p ====
data <- read.table("analysis/006_colorectal_proteins/001_MR_results.txt", header = T, sep = "\t")
data <- subset(data, method == "IVW-MRE")
data <- subset(data, pval < 0.05/1293)

# split into sex and exposure
data_male <- subset(data, group == "Male")
data_colon_male <- subset(data_male, exposure == "colon")
data_distal_male <- subset(data_male, exposure == "distal")
data_proximal_male <- subset(data_male, exposure == "proximal")
data_rectal_male <- subset(data_male, exposure == "rectal")
data_overall_male <- subset(data_male, exposure == "overall")

data_female <- subset(data, group == "Female")
data_colon_female <- subset(data_female, exposure == "colon")
data_distal_female <- subset(data_female, exposure == "distal")
data_proximal_female <- subset(data_female, exposure == "proximal")
data_rectal_female <- subset(data_female, exposure == "rectal")
data_overall_female <- subset(data_female, exposure == "overall")

data_combined <- subset(data, group == "Sex-combined")
data_colon_combined <- subset(data_combined, exposure == "colon")
data_distal_combined <- subset(data_combined, exposure == "distal")
data_proximal_combined <- subset(data_combined, exposure == "proximal")
data_rectal_combined <- subset(data_combined, exposure == "rectal")
data_overall_combined <- subset(data_combined, exposure == "overall")
data_overall_HRC_combined <- subset(data_combined, exposure == "overall-HRC-EUR")
data_overall_early_onset <- subset(data_combined, exposure == "early_onset")

# load directionally associated dataframes 
directionally_associated_colon_male <- read.table("analysis/006_colorectal_proteins/002_colon_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_colon_male <- directionally_associated_colon_male$sequence_id
directionally_associated_distal_male <- read.table("analysis/006_colorectal_proteins/002_distal_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_distal_male <- directionally_associated_distal_male$sequence_id
directionally_associated_rectal_male <- read.table("analysis/006_colorectal_proteins/002_rectal_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_rectal_male <- directionally_associated_rectal_male$sequence_id
directionally_associated_proximal_male <- read.table("analysis/006_colorectal_proteins/002_proximal_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_proximal_male <- directionally_associated_proximal_male$sequence_id
directionally_associated_overall_male <- read.table("analysis/006_colorectal_proteins/002_overall_male_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_overall_male <- directionally_associated_overall_male$sequence_id

directionally_associated_colon_female <- read.table("analysis/006_colorectal_proteins/002_colon_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_colon_female <- directionally_associated_colon_female$sequence_id
directionally_associated_distal_female <- read.table("analysis/006_colorectal_proteins/002_distal_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_distal_female <- directionally_associated_distal_female$sequence_id
directionally_associated_rectal_female <- read.table("analysis/006_colorectal_proteins/002_rectal_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_rectal_female <- directionally_associated_rectal_female$sequence_id
directionally_associated_proximal_female <- read.table("analysis/006_colorectal_proteins/002_proximal_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_proximal_female <- directionally_associated_proximal_female$sequence_id
directionally_associated_overall_female <- read.table("analysis/006_colorectal_proteins/002_overall_female_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_overall_female <- directionally_associated_overall_female$sequence_id

directionally_associated_colon_combined <- read.table("analysis/006_colorectal_proteins/002_colon_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_colon_combined <- directionally_associated_colon_combined$sequence_id
directionally_associated_distal_combined <- read.table("analysis/006_colorectal_proteins/002_distal_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_distal_combined <- directionally_associated_distal_combined$sequence_id
directionally_associated_rectal_combined <- read.table("analysis/006_colorectal_proteins/002_rectal_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_rectal_combined <- directionally_associated_rectal_combined$sequence_id
directionally_associated_proximal_combined <- read.table("analysis/006_colorectal_proteins/002_proximal_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_proximal_combined <- directionally_associated_proximal_combined$sequence_id
directionally_associated_overall_combined <- read.table("analysis/006_colorectal_proteins/002_overall_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_overall_combined <- directionally_associated_overall_combined$sequence_id
directionally_associated_overall_HRC_combined <- read.table("analysis/006_colorectal_proteins/002_overall_HRC_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_overall_HRC_combined <- directionally_associated_overall_HRC_combined$sequence_id
directionally_associated_overall_early_onset <- read.table("analysis/006_colorectal_proteins/002_overall_early_onset_combined_directionally_associated_proteins.txt", header = T, sep = "\t")
directionally_associated_overall_early_onset <- directionally_associated_overall_early_onset$sequence_id

# pull out directionally consistent from p value a data frames
data_colon_male <- data_colon_male[data_colon_male$sequence_id %in% directionally_associated_colon_male, ]
data_distal_male <- data_distal_male[data_distal_male$sequence_id %in% directionally_associated_distal_male, ]
data_proximal_male <- data_proximal_male[data_proximal_male$sequence_id %in% directionally_associated_proximal_male, ]
data_rectal_male <- data_rectal_male[data_rectal_male$sequence_id %in% directionally_associated_rectal_male, ]
data_overall_male <- data_overall_male[data_overall_male$sequence_id %in% directionally_associated_overall_male, ]

data_colon_female <- data_colon_female[data_colon_female$sequence_id %in% directionally_associated_colon_female, ]
data_distal_female <- data_distal_female[data_distal_female$sequence_id %in% directionally_associated_distal_female, ]
data_proximal_female <- data_proximal_female[data_proximal_female$sequence_id %in% directionally_associated_proximal_female, ]
data_rectal_female <- data_rectal_female[data_rectal_female$sequence_id %in% directionally_associated_rectal_female, ]
data_overall_female <- data_overall_female[data_overall_female$sequence_id %in% directionally_associated_overall_female, ]

data_colon_combined <- data_colon_combined[data_colon_combined$sequence_id %in% directionally_associated_colon_combined, ]
data_distal_combined <- data_distal_combined[data_distal_combined$sequence_id %in% directionally_associated_distal_combined, ]
data_proximal_combined <- data_proximal_combined[data_proximal_combined$sequence_id %in% directionally_associated_proximal_combined, ]
data_rectal_combined <- data_rectal_combined[data_rectal_combined$sequence_id %in% directionally_associated_rectal_combined, ]
data_overall_combined <- data_overall_combined[data_overall_combined$sequence_id %in% directionally_associated_overall_combined, ]
data_overall_HRC_combined <- data_overall_HRC_combined[data_overall_HRC_combined$sequence_id %in% directionally_associated_overall_HRC_combined, ]
data_overall_early_onset <- data_overall_early_onset[data_overall_early_onset$sequence_id %in% directionally_associated_overall_early_onset, ]

# make table of N associations
table_p3 <- data.frame(
  exposure = c("overall", "overall_HRC", "early_onset", "colon", "distal", "proximal", "rectal", "unique"),
  male_phenospd_consistent = c(nrow(data_overall_male), 
                                 NA, 
                               NA,
                                 nrow(data_colon_male),
                                 nrow(data_distal_male), 
                                 nrow(data_proximal_male), 
                                 nrow(data_rectal_male), 
                                 length(unique(data_male$sequence_id))),
  
  female_phenospd_consistent = c(nrow(data_overall_female), 
                                   NA, 
                                 NA,
                                   nrow(data_colon_female),
                                   nrow(data_distal_female), 
                                   nrow(data_proximal_female), 
                                   nrow(data_rectal_female), 
                                   length(unique(data_female$sequence_id))),
  
  combined_phenospd_consistent = c(nrow(data_overall_combined), 
                                     nrow(data_overall_HRC_combined), 
                                   nrow(data_overall_early_onset), 
                                   nrow(data_colon_combined),
                                     nrow(data_distal_combined), 
                                     nrow(data_proximal_combined), 
                                     nrow(data_rectal_combined), 
                                     length(unique(data_combined$sequence_id)))
  
)

table <- left_join(table, table_p3)
write.table(table, "analysis/006_colorectal_proteins/table_N_p_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make table of associations ====
data_associations <- bind_rows(
  data_overall_male, data_colon_male, data_distal_male, data_proximal_male, data_rectal_male,
  data_overall_female, data_colon_female, data_distal_female, data_proximal_female, data_rectal_female,
  data_overall_combined, data_overall_HRC_combined, data_overall_early_onset, data_colon_combined, data_distal_combined, data_proximal_combined, data_rectal_combined)

write.table(data_associations, "analysis/006_colorectal_proteins/003_phenospd_directionally_associated_proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# unqiue associations table
data <- read.table("analysis/006_colorectal_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")

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

write.table(table, "analysis/006_colorectal_proteins/table_N_unique_associations.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


