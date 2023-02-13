rm(list=ls())

filenames <- dir("analysis/009_colocalisation/results/", recursive = TRUE, full.names = TRUE, pattern = "001_coloc_results")

coloc_results <- lapply(filenames, fread, header = T)
coloc_results <- bind_rows(coloc_results)

length(unique(coloc_results$exposure))

a <- subset(coloc_results, h4 >= 0.9)

write.table(coloc_results, "analysis/009_colocalisation/results/000_coloc_results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
