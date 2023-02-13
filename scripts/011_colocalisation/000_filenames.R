rm(list=ls())

# make list of files to extract cis_window
filenames <- read.table("data/000_mvmr_data/000_file_names_for_mvmr.txt", header = T, sep = "\t")
filenames$file <- paste0(filenames$file, ".unzipped.cis.txt")

write.table(filenames, "data/000_colocalisation/000_file_names_for_colocalisation.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

