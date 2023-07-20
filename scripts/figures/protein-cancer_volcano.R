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
protein_cancer <- read.table("analysis/011_protein_colorectal_cissnp/001_MR_results.txt", header = T, sep = "\t")
protein_cancer_associations <- read.table("analysis/007_adiposity_proteins_colorectal/protein_cancer.txt", header = T, sep = "\t")
protein_cancer_associations$label <- protein_cancer_associations$protein
plot_data <- left_join(protein_cancer, protein_cancer_associations)
  
# plot data ====
psig <- 0.05/1293
psig_log <- -log10(psig)
table(plot_data$outcome)
plot_data$colour <- ifelse(plot_data$pval < psig, palette_discrete[4], "black")
cancer <- c("overall", "colon", "distal", "proximal", "rectal")
plot_data <- plot_data[plot_data$outcome %in% cancer, ]
plot_data$group <- factor(plot_data$group, levels = c("Sex-combined", "Female", "Male"))
plot_data$outcome <- factor(plot_data$outcome, levels = c("overall", "colon", "proximal", "distal", "rectal"))


# volcano plot ====
tiff(filename = "analysis/figures/protein-crc_volcanoplot.tiff",
     width = 1000, height = 1000, units = "px", res = 100)

ggplot(plot_data, aes(x = b, y = -log10(pval),
                      colour = colour)) +
  geom_point() +
  scale_color_manual(values = c(palette_discrete[4], "black"), guide = FALSE) +
  
  geom_text_repel(aes(label = label), color = "black", 
                  max.iter = 10000, max.overlaps = Inf,
                  min.segment.length = 0,
                  nudge_y = 5, direction = "x",
                  angle = 45) +
  
  facet_grid(group ~ outcome,
             scales = "fixed",
             space = "fixed") +
  theme_cowplot() + 
  theme(legend.position = "bottom") +
  labs(x = "effect estimate", y = "-log10(P-Value)") +
  coord_cartesian(xlim = c(-3,3), ylim = c(0,45))

dev.off()

