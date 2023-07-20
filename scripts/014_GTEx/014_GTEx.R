# gtex download/format - from: https://github.com/tnieuwe/Matrisome_GTEx_analysis
rm(list=ls())

# source ====
library(R.utils)
library(tidyverse)
library(DESeq2)
library(data.table)
library(tibble)
library(dplyr)
library(ggforce)
library(cowplot)

#  data  ====
load("/data/GWAS_data/files/references/GTEx/gtex-gene-counts-v8.rda")

gene <- "GREM1"
df_list <- list()

for (i in 1:length(gene)) {
  
  # filter data for genes of interest
  ind <- which(gtab$gene_name%in% gene[i])
  filt_dat <- data.frame(t(dat[ind,]))
  filt_dat <- rownames_to_column(filt_dat, var = "SAMPID")
  filt_dat <- left_join(filt_dat, stab, by = "SAMPID") 
  colnames(filt_dat)[2] <- "TPM"
  
  # Assign the dataframe to a list element with variable name based on the gene
  assign(paste0(gene[i]), filt_dat, envir = .GlobalEnv)
  df_list[[i]] <- get(paste0(gene[i]))  # Store the dataframe in the list
}

# make sex-specific wide data ====
data_combined <- df_list[[1]]
data_male <- subset(df_list[[1]], SEX == 1)
data_female <- subset(df_list[[1]], SEX == 2)

## combined wide
df <- select(data_combined, SUBJID, SMTSD, TPM)
df_wide <- df %>%
  pivot_wider(names_from = SMTSD, values_from = TPM, id_cols = SUBJID)
df_wide <- as.data.frame(df_wide)
rownames(df_wide) <- df_wide[, 1]  # Set column 1 as row names
data_combined <- df_wide[, -1]  # Remove the first column

## male wide
df <- select(data_male, SUBJID, SMTSD, TPM)
df_wide <- df %>%
  pivot_wider(names_from = SMTSD, values_from = TPM, id_cols = SUBJID)
df_wide <- as.data.frame(df_wide)
rownames(df_wide) <- df_wide[, 1]  # Set column 1 as row names
data_male <- df_wide[, -1]  # Remove the first column

## female wide
df <- select(data_female, SUBJID, SMTSD, TPM)
df_wide <- df %>%
  pivot_wider(names_from = SMTSD, values_from = TPM, id_cols = SUBJID)
df_wide <- as.data.frame(df_wide)
rownames(df_wide) <- df_wide[, 1]  # Set column 1 as row names
data_female <- df_wide[, -1]  # Remove the first column

# run wilcox test ====
## Perform pairwise Wilcoxon rank sum test
reference_col <- "Whole Blood"

## Perform Wilcoxon rank sum test for each column against the reference column

### combined
result <- list()
test_data <- data_combined
for (i in colnames(test_data)) {
  if (i != reference_col) {
    test_result <- wilcox.test(test_data[[i]], test_data[[reference_col]])
    result[[i]] <- test_result$p.value
  }
}
result_combined <- data.frame(
  SMTSD = names(result),
  p_value = unlist(result),
  stringsAsFactors = FALSE
)

### male
result <- list()
test_data <- data_male
for (i in colnames(test_data)) {
  if (i != reference_col) {
    test_result <- wilcox.test(test_data[[i]], test_data[[reference_col]])
    result[[i]] <- test_result$p.value
  }
}
result_male <- data.frame(
  SMTSD = names(result),
  p_value = unlist(result),
  stringsAsFactors = FALSE
)

### female
result <- list()
test_data <- data_female
for (i in colnames(test_data)) {
  if (i != reference_col) {
    test_result <- wilcox.test(test_data[[i]], test_data[[reference_col]])
    result[[i]] <- test_result$p.value
  }
}
result_female <- data.frame(
  SMTSD = names(result),
  p_value = unlist(result),
  stringsAsFactors = FALSE
)

# combine data with stats data ====
result_male$SEX <- 1
result_female$SEX <- 2
data_male <- left_join(df_list[[1]], result_male, by = c("SMTSD", "SEX"))
data_male <- subset(data_male, SEX == 1)
data_female <- left_join(df_list[[1]], result_female, by = c("SMTSD", "SEX"))
data_female <- subset(data_female, SEX == 2)
data_sex <- bind_rows(data_male, data_female)
data_sex <- data_sex %>%
  mutate(SEX = ifelse(SEX == 1, "Male", "Female"))
tissue_male <- unique(data_male$SMTSD)
tissue_female <- unique(data_female$SMTSD)
tissue_combined <- intersect(tissue_male, tissue_female)

data_combined <- left_join(df_list[[1]], result_combined, by = "SMTSD")
data_combined$SEX <- "Sex-combined"
data_combined <- data_combined[data_combined$SMTSD %in% tissue_combined, ]

# plot violin ====
## combined
plot_data <- bind_rows(data_combined,data_sex)
plot_data$SEX <- factor(plot_data$SEX, levels = c("Sex-combined", "Female", "Male"))
plot_data$SMTSD <- fct_relevel(plot_data$SMTSD, "Whole Blood")

ggplot(data = plot_data, 
       aes(x = factor(SMTSD), y = log(TPM),
           colour = SMTS, fill = SMTS)) + 
  geom_violin(width = 1.8, alpha = 0.8) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.8, outlier.alpha = 0) +
  
  labs(xlab = "", ylab = "log(TMP+1)") +
  
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  facet_col(facets = ~SEX, )

