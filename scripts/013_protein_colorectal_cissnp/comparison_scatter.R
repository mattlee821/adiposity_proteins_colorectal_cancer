rm(list=ls())

# environment ====
library(dplyr)
library(data.table)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(tidyr)
library(ggforestplot)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
discrete_palette2 <- wes_palette(colours[3], type = "discrete")

# data ====
a <- read.table("analysis/011_protein_colorectal_cissnp/002_MR_results_shared_with_MVMR.txt", header = T, sep = "\t")
a <- a$ID

mr_all <- fread("analysis/003_proteins_colorectal/001_MR_results.txt")
mr_all$ID <- paste0(mr_all$sequence_id, "_", mr_all$outcome, "_", mr_all$group)
mr_all <- subset(mr_all, method == "IVW-MRE")
mr_all <- mr_all[mr_all$ID %in% a,]

mr_cis <- fread("analysis/011_protein_colorectal_cissnp/001_MR_results.txt")
mr_cis$ID <- paste0(mr_cis$sequence_id, "_", mr_cis$outcome, "_", mr_cis$group)
mr_cis <- mr_cis[mr_cis$ID %in% a,]

data_wide <- left_join(mr_all, mr_cis, by = c("exposure", "group", "outcome"))


# scatter plot ====
plot_data <- data_wide

p1 <- ggplot(data = plot_data, 
             aes(x = b.x, y = b.y, label = exposure)) +
  geom_point(aes(colour = outcome, shape = group), size = 5)  +
  
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
  xlab("beta cis- and trans-SNP (IVW-MRE)") + 
  ylab("beta cis-SNP (Wald ratio)") 

p1

pdf("analysis/011_protein_colorectal_cissnp/figures/shared.pdf",
    width = 10, height = 8 , pointsize = 10)
p1
dev.off()
