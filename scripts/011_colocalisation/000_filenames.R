rm(list=ls())

# make list of files to extract cis_window
filenames <- read.table("data/000_mvmr_data/000_file_names_for_mvmr.txt", header = T, sep = "\t")
filenames$file <- paste0(filenames$file, ".unzipped.cis.txt")

filenames <- as.data.frame(list.files(path = "/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb/", full.names = TRUE))
colnames(filenames) <- "file"
filenames$file <- gsub("/data/protein_GWAS_ferkingstad_EU_2021/files/cis_snps_1mb//", "", filenames$file)

write.table(filenames, "data/000_colocalisation/000_file_names_for_colocalisation.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

