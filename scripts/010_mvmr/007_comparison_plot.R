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
colours <- names(wes_palettes)
discrete_palette <- wes_palette(colours[8], type = "discrete")
source("scripts/my_forestplot.R")

# mvmr data ====
mvmr <- read.table("analysis/008_mvmr/mvmr_results.txt", header = T, sep = "\t")
mvmr$ID_adiposity_cancer <- paste0(mvmr$adiposity, "_", mvmr$cancer, "_", mvmr$group)
ID_adiposity_cancer_mvmr <- unique(mvmr$ID_adiposity_cancer)
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
mr <- subset(mr, method == "IVW-MRE")
mr$ID_adiposity_cancer <- paste0(mr$exposure, "_", mr$outcome, "_", mr$group)
mr <- mr[mr$ID_adiposity_cancer %in% ID_adiposity_cancer_mvmr, ]
mr <- mr[,c("ID_adiposity_cancer", "b", "se", "pval", "exposure", "exposure", "outcome", "group")]
mr$method <- "UMR"
colnames(mr) <- c("ID_adiposity_cancer", "b", "se", "p", "exposure", "adiposity", "cancer", "group", "method")
mr$protein <- "Univariable MR"
mr$ID1 <- mr$exposure
mr$order <- as.factor(0)
mr$x_axis <- mr$exposure

# join ====
data <- bind_rows(mr,mvmr)
data$cancer <- as.factor(data$cancer)
data$cancer <- factor(data$cancer, levels = c("overall", "overall-HRC-EUR", "colon", "distal", "proximal", "rectal"))
data$group <- factor(data$group, levels = c("Sex-combined", "Male", "Female"))
data$method <- factor(data$method, levels = c("UMR", "MVMR"))
data <- data[order(data$order),]

# join proportion mediated ====
proportion_mediated <- read.table("analysis/008_mvmr/proportion_mediated.txt", header = T, sep = "\t")
data <- left_join(data, proportion_mediated)


# data for plot ====
data_bmi <- subset(data, adiposity == "BMI")
data_whr <- subset(data, adiposity == "WHR")

# plot variables====
psignif <- 1
ci <- 0.95

# bmi plot ====
## overall ====
plot_data <- data_bmi
plot_data <- subset(plot_data, cancer == "overall")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Overall colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/BMI_overall.pdf",
    width = 13, height = 8 , pointsize = 10)
p1
dev.off()

## colon ====
plot_data <- data_bmi
plot_data <- subset(plot_data, cancer == "colon")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Colon colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/BMI_colon.pdf",
    width = 13, height = 9, pointsize = 10)
p1
dev.off()

## proximal ====
plot_data <- data_bmi
plot_data <- subset(plot_data, cancer == "proximal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.038)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Proximal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/BMI_proximal.pdf",
    width = 13, height = 24, pointsize = 10)
p1
dev.off()

## distal ====
plot_data <- data_bmi
plot_data <- subset(plot_data, cancer == "distal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Distal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/BMI_distal.pdf",
    width = 13, height = 9, pointsize = 10)
p1
dev.off()

## rectal ====
plot_data <- data_bmi
plot_data <- subset(plot_data, cancer == "rectal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Rectal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/BMI_rectal.pdf",
    width = 13, height = 10, pointsize = 10)
p1
dev.off()


# WHR plot ====
## overall ====
plot_data <- data_whr
plot_data <- subset(plot_data, cancer == "overall")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Overall colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/WHR_overall.pdf",
    width = 13, height = 5, pointsize = 10)
p1
dev.off()

## colon ====
plot_data <- data_whr
plot_data <- subset(plot_data, cancer == "colon")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Colon colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/WHR_colon.pdf",
    width = 13, height = 4 , pointsize = 10)
p1
dev.off()

## proximal ====
plot_data <- data_whr
plot_data <- subset(plot_data, cancer == "proximal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Proximal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/WHR_proximal.pdf",
    width = 13, height = 8, pointsize = 10)
p1
dev.off()

## distal ====
plot_data <- data_whr
plot_data <- subset(plot_data, cancer == "distal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Distal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/WHR_distal.pdf",
    width = 13, height = 3, pointsize = 10)
p1
dev.off()

## rectal ====
plot_data <- data_whr
plot_data <- subset(plot_data, cancer == "rectal")
plot_data <- plot_data[plot_data$ID1 %in% adiposity_labels, ]
plot_data$label <- paste0(round(exp(plot_data$b),2), "; ", round(exp(plot_data$b - (1.96 * plot_data$se)),2), " - ", round(exp(plot_data$b + (1.96 * plot_data$se)),2))

x_min <- min(exp(plot_data$b - (1.96 * plot_data$se)))
x_max <- max(exp(plot_data$b + (1.96 * plot_data$se)))

p1 <- my_forestplot(df = plot_data,
                    name = x_axis,
                    estimate = b,
                    pvalue = p,
                    psignif = psignif,
                    ci = ci,
                    se = se,
                    colour = method,
                    logodds = T) +
  xlab(label = "Odds ratio and 95% confidence interval") +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(1, x_max*1.035)) +
  scale_x_continuous() +
  scale_color_manual(values = c(discrete_palette[3], discrete_palette[1])) +
  ggtitle("Rectal colorectal cancer") +
  ggforce::facet_col(facets = ~group,
                     scales = "free_y",
                     space = "free") +
  aes(label = plot_data$label) +
  geom_text(x = x_max*1.01, vjust = 0.5, hjust = 0)

pdf("analysis/008_mvmr/figures/WHR_rectal.pdf",
    width = 13, height = 5, pointsize = 10)
p1
dev.off()

