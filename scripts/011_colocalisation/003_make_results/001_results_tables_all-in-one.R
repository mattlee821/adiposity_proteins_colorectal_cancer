directories <- list.dirs("analysis/009_colocalisation/results/", recursive = F)
directories <- gsub("analysis/009_colocalisation/results//", "", directories)

for (i in 1:length(directories)){
  outcome <- directories[i]
  
  path <- "analysis/009_colocalisation/results/"
  filenames <- list.files(path = paste0(path, outcome, "/"), pattern = ".RData", full.names = TRUE)
  
  exposure <- gsub(".RData", "", filenames)
  exposure <- gsub(paste0("analysis/009_colocalisation/results/", outcome, "//"), "", exposure)
  exposure <- gsub(paste0("_",outcome), "", exposure)
  
  table_master <- data.frame() # make empty dataframe for final results
  
  for (j in 1:length(filenames)){
    coloc_results <- readRDS(filenames[j])
    
    label_exposure <- exposure[j]
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
  
}

filenames <- list.files(path = path, pattern = "002_coloc_results.txt", recursive = T, full.names = TRUE)
coloc_results <- lapply(filenames, fread)
coloc_results <- bind_rows(coloc_results)
coloc_results <- subset(coloc_results, outcome != "overall_HRC")

coloc_results <- coloc_results %>%
  mutate(cancer = case_when(
    grepl('colon', outcome) ~ 'colon',
    grepl('distal', outcome) ~ 'distal',
    grepl('proximal', outcome) ~ 'proximal',
    grepl('rectal', outcome) ~ 'rectal',
    grepl('stratified_female', outcome) ~ 'overall',
    grepl('stratified_male', outcome) ~ 'overall',
    grepl('overall', outcome) ~ 'overall',
    TRUE ~ 'unknown'
  ))

coloc_results <- coloc_results %>%
  mutate(group = case_when(
    grepl('female', outcome) ~ 'Female',
    grepl('male', outcome) ~ 'Male',
    grepl('overall', outcome) ~ 'Sex-combined',
    grepl('stratified_colon', outcome) ~ 'Sex-combined',
    grepl('stratified_distal', outcome) ~ 'Sex-combined',
    grepl('stratified_proximal', outcome) ~ 'Sex-combined',
    grepl('stratified_rectal', outcome) ~ 'Sex-combined',
    TRUE ~ 'unknown'
  ))

coloc_results <- select(coloc_results, exposure, cancer, group, id, nsnps, h0, h1, h2, h3, h4)
colnames(coloc_results)[2] <- "outcome"
coloc_results$id <- paste0(coloc_results$exposure, "_", coloc_results$outcome, "_", coloc_results$group)
write.table(coloc_results, paste0(path, "/002_coloc_results.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

a <- subset(coloc_results, h4 >= 0.8)
table(a$outcome, a$group)
b <- unique(a$exposure)

coloc_results <- read.table("analysis/009_colocalisation/results/002_coloc_results.txt", header = T, sep = "\t")
