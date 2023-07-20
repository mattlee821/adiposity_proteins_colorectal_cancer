rm(list=ls())

# environment ====
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gghighlight)
library(functions)
library(wesanderson)
library(dplyr)
palette_discrete <- palette()

# data ====
adiposity_protein <- read.table("analysis/002_adiposity_proteins/001_MR_results.txt", header = T, sep = "\t")
adiposity_protein <- subset(adiposity_protein, method == "IVW-MRE")
adiposity_protein_associatons <- read.table("analysis/002_adiposity_proteins/003_phenospd_directionally_associated_proteins.txt", header = T, sep = "\t")
adiposity_protein_associatons <- select(adiposity_protein_associatons, exposure, protein, group)
adiposity_protein_associatons$label <- adiposity_protein_associatons$protein
plot_data <- left_join(adiposity_protein, adiposity_protein_associatons)
table(plot_data$exposure)
adiposity <- c("BMI", "WHR")
plot_data <- plot_data[plot_data$exposure %in% adiposity, ]

protein_cancer_associations <- read.table("analysis/007_adiposity_proteins_colorectal/protein_cancer.txt", header = T, sep = "\t")
protein_cancer_associations <- select(protein_cancer_associations, sequence_id, group, protein)
colnames(protein_cancer_associations)[3] <- "label2"
protein_cancer_associations <- unique(protein_cancer_associations)
plot_data <- left_join(plot_data, protein_cancer_associations)

# plot data ====
psig <- 0.05/1293
psig_log <- -log10(psig)
table(plot_data$exposure)

plot_data$colour <- ifelse(plot_data$pval > 0.05 & is.na(plot_data$label), "black",
                        ifelse(plot_data$pval > 0.05 & !is.na(plot_data$label), "black",
                               ifelse(plot_data$pval < 0.05 & is.na(plot_data$label), "black", # if just sig
                                      palette_discrete[6]))) # if sig and labeled

plot_data$group <- factor(plot_data$group, levels = c("Sex-combined", "Female", "Male"))


# volcano plot ====
tiff(filename = "analysis/figures/adiposity-protein_volcanoplot.tiff",
     width = 1000, height = 1000, units = "px", res = 100)

ggplot(plot_data, aes(x = b, y = -log10(pval),
                      colour = colour)) +
  geom_point() +
  scale_color_manual(values = c(palette_discrete[6], "black"),
                     labels = c("adiposity-protein associations", ""),
                     drop = T) +
  guides(color = guide_legend(override.aes = list(shape = c(19, NA)))) +
  
  geom_text_repel(aes(label = label2), color = "black", 
                  max.iter = 10000, max.overlaps = Inf,
                  min.segment.length = 0,
                  ylim = 50, nudge_y = 10, direction = "x",
                  angle = 45) +

  facet_grid(group ~ exposure,
             scales = "fixed",
             space = "fixed") +
  theme_cowplot() + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "effect estimate", y = "-log10(P-Value)") 

dev.off()
