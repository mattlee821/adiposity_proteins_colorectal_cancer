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
library(ggpubr)
library(ggrepel)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
source("scripts/my_forestplot.R")
source("scripts/my_mr_scatter_plot.R")

# shared ====
data <- read.table("analysis/010_UKB_proteins_colorectal/001_MR_results.txt", header = T, sep = "\t")
data <- subset(data, method == "IVW-MRE")
proteins_UKB <- unique(data$exposure)
data <- read.table("analysis/008_mvmr/mvmr_results.txt", heade = T, sep = "\t")
data <- subset(data, mediator != "BMI")
data <- subset(data, mediator != "WHR")
proteins_mvmr <- unique(data$protein)
proteins_shared <- intersect(proteins_UKB, proteins_mvmr)
data <- data[data$protein %in% proteins_shared,]
sequenceID_mvmr <- unique(data$mediator)

# mvmr ====
data <- read.table("analysis/008_mvmr/mvmr_results.txt", heade = T, sep = "\t")
data <- data[data$mediator %in% sequenceID_mvmr,]
ID <- unique(paste0(data$protein, "_", data$cancer, "_", data$group))

# extract shared analyses ====
## ferkingstad ====
data_ferkingstad <- read.table("analysis/003_proteins_colorectal/001_MR_results.txt", header = T, sep = "\t")
data_ferkingstad <- data_ferkingstad[data_ferkingstad$sequence_id %in% sequenceID_mvmr,]
data_ferkingstad <- subset(data_ferkingstad, method == "IVW-MRE")
data_ferkingstad$ID <- paste0(data_ferkingstad$protein, "_", data_ferkingstad$outcome, "_", data_ferkingstad$group)
data_ferkingstad <- data_ferkingstad[data_ferkingstad$ID %in% ID,]
data_ferkingstad$analysis <- "Ferkingstad"

## UKB ====
data_UKB <- read.table("analysis/010_UKB_proteins_colorectal/001_MR_results.txt", header = T, sep = "\t")
data_UKB$ID <- paste0(data_UKB$exposure, "_", data_UKB$outcome, "_", data_UKB$group)
data_UKB <- data_UKB[data_UKB$ID %in% ID,]
data_UKB <- subset(data_UKB, method == "IVW-MRE")
data_UKB$analysis <- "Sun"

# combine ====
data_ferkingstad <- select(data_ferkingstad, protein, outcome, group, nsnp, method, b, se, pval, OR, lower_ci, upper_ci, analysis)
colnames(data_ferkingstad)[1] <- "exposure"
data_UKB <- select(data_UKB, exposure, outcome, group, nsnp, method, b, se, pval, OR, lower_ci, upper_ci, analysis)
data <- rbind(data_ferkingstad, data_UKB)

data_long <- rbind(data_UKB, data_ferkingstad)

data_wide <- left_join(data_UKB, data_ferkingstad, by = c("exposure", "group", "outcome"))
data_wide$outcome <- factor(data_wide$outcome, levels = c("overall",
                                                          "colon", "distal", "proximal", "rectal"))

# scatter plot ====
plot_data <- data_wide

p1 <- ggplot(data = plot_data, 
       aes(x = b.x, y = b.y, label = exposure)) +
  geom_point(aes(colour = outcome, shape = group), size = 5)  +
  
  geom_text_repel(point.padding = 10) +

  geom_errorbar(aes(ymin = b.y - se.y, 
                    ymax = b.y + se.y), 
                colour = "grey", width = 0) +
  geom_errorbarh(aes(xmin = b.x - se.x, 
                     xmax = b.x + se.x),
                 colour = "grey", height = 0) + 
  
  guides(colour = guide_legend(override.aes = list(size = 2),
                               title = "",
                               label.hjust = 0,
                               label.vjust = 0.5)) +
  
  guides(shape = guide_legend(override.aes = list(size = 2),
                               title = "",
                               label.hjust = 0,
                               label.vjust = 0.5)) +
  
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +

  theme_cowplot() +
  xlab("beta Ferkingstad et al., proteins > cancer") + 
  ylab("beta Sun et al., proteins > cancer") 

p1

pdf("analysis/010_UKB_proteins_colorectal/figures/shared.pdf",
    width = 10, height = 8 , pointsize = 10)
p1
dev.off()

