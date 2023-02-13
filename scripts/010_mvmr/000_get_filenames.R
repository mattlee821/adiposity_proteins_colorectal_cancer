rm(list=ls())
library(dplyr)

a <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/003_adiposity_protein_cancer_pos.txt", header = T, sep = "\t")
a <- paste0(a$sequence_id, "_", a$gene, "_", a$protein, ".txt.gz.annotated.gz.exclusions.gz.alleles.gz")
id <- as.data.frame(unique(a))
colnames(id) <- "file"

a <- read.table("analysis/007_adiposity_proteins_colorectal/phenospd/003_adiposity_protein_cancer_neg.txt", header = T, sep = "\t")
a <- paste0(a$sequence_id, "_", a$gene, "_", a$protein, ".txt.gz.annotated.gz.exclusions.gz.alleles.gz")
id1 <- as.data.frame(unique(a))
colnames(id1) <- "file"

id <- bind_rows(id,id1)

write.table(id, "data/000_mvmr_data/000_file_names_for_mvmr.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
