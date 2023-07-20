rm(list=ls())

outcome <- "joint_colon_female"
path <- "analysis/009_colocalisation/results/"
filenames <- list.files(path = paste0(path, outcome, "/"), pattern = ".RData", full.names = TRUE)

exposure <- gsub(".RData", "", filenames)
exposure <- gsub("analysis/009_colocalisation/results/joint_colon_female//", "", exposure)
exposure <- gsub(paste0("_",outcome), "", exposure)

table_master <- data.frame() # make empty dataframe for final results

for (i in 1:length(filenames)){
  coloc_results <- readRDS(filenames[i])
  
  label_exposure <- exposure[i]
  label_outcome <- outcome
  label <- paste0(label_exposure, "_", label_outcome)
  
  # make table ====
  table <- data.frame(
    exposure = label_exposure,
    outcome = label_outcome,
    id = label,
    nsnps = coloc_results["summary"][[1]][1],
    h0 = coloc_results["summary"][[1]][2],
    h1 = coloc_results["summary"][[1]][3],
    h2 = coloc_results["summary"][[1]][4],
    h3 = coloc_results["summary"][[1]][5],
    h4 = coloc_results["summary"][[1]][6])
  
  table_master <- rbind(table_master, table)
  
}

write.table(table_master, paste0(path, outcome, "/002_coloc_results.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

