rm(list=ls())
library(dplyr)

a <- read.table("analysis/007_adiposity_proteins_colorectal/adiposity_protein_cancer.txt", header = T, sep = "\t")
a <- paste0(a$sequence_id, "_", a$gene, "_", a$protein, ".txt.gz.annotated.gz.exclusions.gz.alleles.gz")
id <- as.data.frame(unique(a))
colnames(id) <- "file"

write.table(id, "data/000_mvmr_data/000_file_names_for_mvmr.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
