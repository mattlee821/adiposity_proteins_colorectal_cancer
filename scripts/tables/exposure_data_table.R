rm(list=ls())

# adiposity ====
data <- fread("data/exposure_data/000_clumped_adiposity.txt", header = T)
head(data)
data$group[data$id.exposure == "BMI_combined"] <- "Sex-combined"
data$group[data$id.exposure == "WHR_combined"] <- "Sex-combined"
data$group[data$id.exposure == "BMI_male"] <- "Male"
data$group[data$id.exposure == "WHR_male"] <- "Male"
data$group[data$id.exposure == "BMI_female"] <- "Female"
data$group[data$id.exposure == "WHR_female"] <- "Female"
data$dataset <- "adiposity"
adiposity <- select(data, dataset, id.exposure, exposure, group, CHR,POS,SNP,effect_allele.exposure,other_allele.exposure, eaf.exposure,beta.exposure, se.exposure, pval.exposure, samplesize.exposure)

# cancer ====
data <- fread("data/exposure_data/000_clumped_crc.txt", header = T)
head(data)
data$group[data$id.exposure == "overall_combined"] <- "Sex-combined"
data$group[data$id.exposure == "overall_combined_HRC"] <- "Sex-combined"
data$group[data$id.exposure == "colon_combined"] <- "Sex-combined"
data$group[data$id.exposure == "distal_combined"] <- "Sex-combined"
data$group[data$id.exposure == "proximal_combined"] <- "Sex-combined"
data$group[data$id.exposure == "rectal_combined"] <- "Sex-combined"

data$group[data$id.exposure == "overall_male"] <- "Male"
data$group[data$id.exposure == "colon_male"] <- "Male"
data$group[data$id.exposure == "distal_male"] <- "Male"
data$group[data$id.exposure == "proximal_male"] <- "Male"
data$group[data$id.exposure == "rectal_male"] <- "Male"

data$group[data$id.exposure == "overall_female"] <- "Female"
data$group[data$id.exposure == "colon_female"] <- "Female"
data$group[data$id.exposure == "distal_female"] <- "Female"
data$group[data$id.exposure == "proximal_female"] <- "Female"
data$group[data$id.exposure == "rectal_female"] <- "Female"


data$samplesize.exposure[data$id.exposure == "overall_combined"] <- 58131+67347
data$samplesize.exposure[data$id.exposure == "overall_combined_HRC"] <- 55168+65160
data$samplesize.exposure[data$id.exposure == "colon_combined"] <- 32002+64159
data$samplesize.exposure[data$id.exposure == "distal_combined"] <- 14376+64159
data$samplesize.exposure[data$id.exposure == "proximal_combined"] <- 15706+64159
data$samplesize.exposure[data$id.exposure == "rectal_combined"] <- 16212+64159

data$samplesize.exposure[data$id.exposure == "overall_male"] <- 31288+34527
data$samplesize.exposure[data$id.exposure == "colon_male"] <- 15897+34527
data$samplesize.exposure[data$id.exposure == "distal_male"] <- 7687+34527
data$samplesize.exposure[data$id.exposure == "proximal_male"] <- 7186+34527
data$samplesize.exposure[data$id.exposure == "rectal_male"] <- 9671+34527

data$samplesize.exposure[data$id.exposure == "overall_female"] <- 26843+32820
data$samplesize.exposure[data$id.exposure == "colon_female"] <- 16105+32820
data$samplesize.exposure[data$id.exposure == "distal_female"] <- 6689+32820
data$samplesize.exposure[data$id.exposure == "proximal_female"] <- 8520+32820
data$samplesize.exposure[data$id.exposure == "rectal_female"] <- 6541+32820

data$dataset <- "cancer"
colnames(data)[3] <- "POS"
colnames(data)[4] <- "CHR"
cancer <- select(data, dataset, id.exposure, exposure, group, CHR,POS,SNP,effect_allele.exposure,other_allele.exposure, eaf.exposure,beta.exposure, se.exposure, pval.exposure, samplesize.exposure)

# proteins decode ====
data <- fread("data/000_protein_data/exposure_data/000_clumped.txt", header = T)
head(data)
data$CHR <- gsub("chr", "", data$CHR)
data$id.exposure <- paste0(data$sequence_id, "_", data$id.exposure)
data$exposure <- data$id.exposure
data$dataset <- "proteins_Ferkingstad"
data$group <- "Sex-combined"
proteins_decode <- select(data, dataset, id.exposure, exposure, group, CHR,POS,SNP,effect_allele.exposure,other_allele.exposure, eaf.exposure,beta.exposure, se.exposure, pval.exposure, samplesize.exposure)

# protein ukb ====
data <- fread("analysis/010_UKB_proteins_colorectal/exposure_data.txt", stringsAsFactors = T)
head(data)
data$dataset <- "proteins_Sun"
data$group <- "Sex-combined"
data$samplesize.exposure <- 35571
colnames(data) <- c("SNP", "CHR", "POS", "effect_allele.exposure", "other_allele.exposure", "exposure", "eaf.exposure", "beta.exposure", "se.exposure", "log10p.exposure","pval.exposure", "id.exposure","f_stats", "dataset","group", "samplesize.exposure")
proteins_ukb <- select(data, dataset, id.exposure, exposure, group, CHR,POS,SNP,effect_allele.exposure,other_allele.exposure, eaf.exposure,beta.exposure, se.exposure, pval.exposure, samplesize.exposure)

# proteins decode cis-SNP ====
data <- fread("data/000_protein_data/exposure_data/cis_snps/cis_snps.txt", header = T)
ncols <- max(stringr::str_count(data$phenotype, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = phenotype,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_ID <- paste0(data$col1, "_", data$col2)
data$id.exposure <- sub(".*?_", "", data$phenotype) # remove first part of sequence_id from column
data$id.exposure <- sub(".*?_", "", data$id.exposure) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.exposure) # remove genename
data$gene <- gsub("\\_.*","",data$id.exposure)
data <- select(data, -col1,-col2,-col3,-col4,-col5,-col6,-col7,-col8)
colnames(data) <- c("CHR", "POS", "ID", "SNP", "A1", "A2", "beta.exposure","pval.exposure","log10p.exposure","se.exposure","samplesize.exposure","ImpMAF","id.exposure",
                    "effect_allele.exposure", "other_allele.exposure", "eaf.exposure",  "sequence_ID", "exposure", "protein","gene")
data$dataset <- "proteins_Ferkingstad_cis"
data$group <- "Sex-combined"
data$CHR <- gsub("chr", "", data$CHR)
data <- select(data, dataset, id.exposure, exposure, group, CHR,POS,SNP,effect_allele.exposure,other_allele.exposure, eaf.exposure,beta.exposure, se.exposure, pval.exposure, samplesize.exposure)
proteins_cis <- data[ , .SD[which.min(pval.exposure)], by = id.exposure]

# master ====
data <- rbind(adiposity, cancer, proteins_decode, proteins_cis, proteins_ukb)
data$f_statistic <- (data$beta.exposure / data$se.exposure)^2 
data <- subset(data, dataset != "cancer")
write.table(data, "analysis/tables/exposure_data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

