rm(list=ls())

# source ====
library(ggplot2)
library(dplyr)
library(knitr)
library(patchwork)
library(tidyr)
library(purrr)
devtools::install_github("mattlee821/ggforestplot")
library(ggforestplot)
library(wesanderson)
library(rlang)
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")

# data ====
data <- read.table("analysis/001_adiposity_colorectal/001_MR_results.txt", header = T, sep = "\t", stringsAsFactors=T)
data <- subset(data, outcome != "overall-HRC-EUR")
# plot_data ====
plot_data <- data 
plot_data$exposure <- factor(plot_data$exposure, levels = c("BMI", "WHR"))
plot_data$outcome <- factor(plot_data$outcome, levels = c(
  "overall", "colon", "proximal", "distal", "rectal"))
plot_data$group <- factor(plot_data$group, levels = c("Sex-combined", "Male", "Female"))
plot_data <- droplevels(plot_data)

plot_data$order[plot_data$outcome == "overall"] <- 1
plot_data$order[plot_data$outcome == "colon"] <- 2
plot_data$order[plot_data$outcome == "distal"] <- 4
plot_data$order[plot_data$outcome == "proximal"] <- 3
plot_data$order[plot_data$outcome == "rectal"] <- 5

xmin <- min(plot_data$lower_ci)
xmax <- max(plot_data$upper_ci)
psignif <- 1
ci <- 0.95

plot_data <- plot_data[order(plot_data$order),]

pdf("analysis/001_adiposity_colorectal/figures/forestplot_all.pdf",
    width = 10, height = 14, pointsize = 10)
forestplot(df = plot_data,
              name = outcome,
              estimate = b,
              pvalue = pval,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = exposure,
              shape = method,
              logodds = T, 
              space = 0.9) +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1], discrete_palette[5])) +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank()) 
dev.off()

plot_data <- subset(plot_data, method == "IVW-MRE")
pdf("analysis/001_adiposity_colorectal/figures/forestplot_main.pdf",
    width = 10, height = 10, pointsize = 10)
forestplot(df = plot_data,
              name = outcome,
              estimate = b,
              pvalue = pval,
              psignif = psignif,
              ci = ci,
              se = se,
              colour = exposure,
              logodds = T, 
              space = 0.9) +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1], discrete_palette[5])) +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()

