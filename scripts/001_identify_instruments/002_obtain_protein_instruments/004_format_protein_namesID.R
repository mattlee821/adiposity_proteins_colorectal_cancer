library(tidyr)
library(data.table)

# clumped data ====
# data ====
data <- read.table("data/exposure_data/000_clumped_proteins.txt", header = T, sep = "\t")

# sequence_id column ====
ncols <- max(stringr::str_count(data$id.exposure, "_")) + 1
colmn <- paste0("col", 1:ncols)
data <-
  tidyr::separate(
    data = data,
    col = id.exposure,
    sep = "_",
    into = colmn,
    remove = FALSE)
data$sequence_id <- paste0(data$col1, "_", data$col2)

data <- data[,c("CHR","POS","SNPID","SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","pval.exposure",
                "minus_log10_pval","se.exposure","samplesize.exposure","eaf.exposure","sequence_id","exposure","id.exposure")]

# protein_name column ====
data$id.exposure <- sub(".*?_", "", data$id.exposure) # remove first part of sequence_id from column
data$id.exposure <- sub(".*?_", "", data$id.exposure) # remove second part of sequence_id from column
data$protein <- sub(".*?_", "", data$id.exposure) # remove genename

# gene name column 
data$gene <- gsub("\\_.*","",data$id.exposure)

# save ====
write.table(data, "data/000_protein_data/exposure_data/000_clumped.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
