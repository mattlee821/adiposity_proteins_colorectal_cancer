rm(list=ls())

# source ====
library(ggplot2)
library(dplyr)
library(knitr)
library(patchwork)
library(tidyr)
library(purrr)
library(ggforestplot)
library(wesanderson)
library(EpiViz)
library(rlang)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
source("scripts/my_forestplot.R")

# adiposity_protein_cancer ====
adiposity_protein_cancer <- read.table("analysis/007_adiposity_proteins_colorectal/adiposity_protein_cancer.txt", header = T, sep = "\t")
adiposity_protein <- paste0(adiposity_protein_cancer$adiposity, "_", adiposity_protein_cancer$outcome, "_", adiposity_protein_cancer$group)
protein_cancer <- paste0(adiposity_protein_cancer$sequence_id, "_", adiposity_protein_cancer$outcome, "_", adiposity_protein_cancer$group)
adiposity_cancer <- paste0(adiposity_protein_cancer$adiposity, "_", adiposity_protein_cancer$outcome, "_", adiposity_protein_cancer$group)
adiposity_protein_cancer <- paste0(adiposity_protein_cancer$adiposity, "_",adiposity_protein_cancer$sequence_id, "_", adiposity_protein_cancer$outcome, "_", adiposity_protein_cancer$group)

# mvmr data ====
mvmr <- read.table("analysis/008_mvmr/mvmr_results.txt", header = T, sep = "\t")
mvmr <- mvmr[mvmr$ID %in% adiposity_protein_cancer, ]
mvmr <- mvmr[,c("ID", "b", "se", "p", "exposure", "adiposity", "protein", "cancer", "group")]
mvmr$method <- "MVMR"
mvmr$x_axis <- mvmr$protein

## subset for mediators
adiposity_labels <- c("BMI", "WHR", "WHRadjBMI")
mvmr_adiposity <- mvmr[mvmr$exposure %in% adiposity_labels, ]
mvmr_adiposity$ID1 <- mvmr_adiposity$exposure
mvmr_protein <- mvmr[!mvmr$exposure %in% adiposity_labels, ]
mvmr_protein$ID1 <- "protein"
mvmr <- bind_rows(mvmr_adiposity, mvmr_protein)
mvmr$order <- factor(mvmr$protein, labels = 1:nlevels( as.factor(mvmr$protein)))

# 2SMR data ====
mr <- read.table("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T, sep = "\t")
mr$ID <- paste0(mr$exposure, "_", mr$outcome, "_", mr$group)
mr <- mr[mr$ID %in% adiposity_cancer, ]
mr <- subset(mr, method == "IVW-MRE")
mr <- mr[,c("ID", "b", "se", "pval", "exposure", "exposure", "outcome", "group")]
mr$method <- "UVMR"
colnames(mr) <- c("ID", "b", "se", "p", "exposure", "adiposity", "cancer", "group", "method")
mr$protein <- "Univariable MR"
mr$ID1 <- mr$exposure
mr$order <- as.factor(0)
mr$x_axis <- mr$exposure

# join ====
data <- bind_rows(mr,mvmr)
data$cancer <- as.factor(data$cancer)
data$cancer <- factor(data$cancer, levels = c("overall", "overall-HRC-EUR", "colon", "distal", "proximal", "rectal"))
data$group <- factor(data$group, levels = c("Sex-combined", "Male", "Female"))
data$method <- factor(data$method, levels = c("UVMR", "MVMR"))
data <- data[order(data$order),]
write.table(data, "analysis/008_mvmr/MVMR_UVMR_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# join proportion mediated ====
proportion_mediated <- read.table("analysis/008_mvmr/proportion_mediated.txt", header = T, sep = "\t")
data <- left_join(data, proportion_mediated)

# plot variables====
psignif <- 1
ci <- 0.95

# plot ====
plot_data <- subset(data, exposure == "BMI")
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))
plot_data$x_axis <- factor(plot_data$x_axis, levels = c("GREM1", "BMI"))
plot_data$cancer <- factor(plot_data$cancer, levels = c("overall", "colon", "distal", "proximal", "rectal"))

levels(plot_data$cancer)

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = cancer,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = x_axis,
                    logodds = T) +
  
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  
  labs(x = "Odds ratio and 95% confidence interval", colour = "") +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[1], discrete_palette[3]), 
                     labels = c("MVMR BMI adjusted for GREM1", "UVMR BMI")) +
  scale_y_discrete(limits = rev(c("overall", "colon", "distal", "proximal", "rectal"))) +
  
  theme(legend.position = "bottom") 

p1 

tiff("analysis/008_mvmr/figures/adiposity_protein_cancer_MVMR.tiff", width = 1000, height = 1000, units = "px", res = 100)
p1
dev.off()
