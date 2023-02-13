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
library(rlang)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
source("scripts/my_forestplot.R")

# data ====
data <- read.table("analysis/004_colorectal_adiposity/001_MR_results.txt", header = T, sep = "\t", stringsAsFactors=T)

# plot_data ====
plot_data <- data 
plot_data$outcome <- factor(plot_data$outcome, levels = c("BMI", "WHR", "WHRadjBMI"))
plot_data$exposure <- factor(plot_data$exposure, levels = c(
  "overall", "overall_HRC", "colon", "proximal", "distal", "rectal"))
plot_data$group <- factor(plot_data$group, levels = c("Sex-combined", "Male", "Female"))
plot_data <- droplevels(plot_data)

plot_data$order[plot_data$exposure == "overall"] <- 1
plot_data$order[plot_data$exposure == "overall_HRC"] <- 2
plot_data$order[plot_data$exposure == "colon"] <- 3
plot_data$order[plot_data$exposure == "distal"] <- 4
plot_data$order[plot_data$exposure == "proximal"] <- 5
plot_data$order[plot_data$exposure == "rectal"] <- 6

xmin <- min(plot_data$lower_ci)
xmax <- max(plot_data$upper_ci)
psignif <- 0.05
ci <- 0.95

plot_data <- plot_data[order(plot_data$order),]

pdf("analysis/004_colorectal_adiposity/figures/forestplot_all.pdf",
    width = 10, height = 14, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = pval,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = outcome,
              shape = method,
              logodds = F, 
              space = 0.9) +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1], discrete_palette[5])) +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

plot_data <- subset(plot_data, method == "IVW-MRE")
pdf("analysis/004_colorectal_adiposity/figures/forestplot_main.pdf",
    width = 10, height = 10, pointsize = 10)
my_forestplot(df = plot_data,
              name = exposure,
              estimate = b,
              pvalue = pval,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = outcome,
              logodds = F, 
              space = 0.9) +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1], discrete_palette[5])) +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()
