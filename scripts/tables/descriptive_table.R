library(data.table)
library(dplyr)

# data tables ====
table_adiposity <- data.frame(
  data_set = c(rep("adiposity", 6)), 
  trait = c("BMI", "BMI", "BMI", "WHR", "WHR", "WHR"),
  sex = c("combined", "female", "male", "combined", "female", "male"),
  population = c(rep("European", 6)),
  N = c("806834", "434794", "374756", "697734", "381152", "316772"),
  N_control = c(rep(NA, 6)),
  pvalue_threshold = c(rep(5E-9, 6)),
  LD_r2 = c(rep(0.001, 6)),
  LD_window = c(rep(10000, 6)),
  SNPs = c(457,245,204, 283,206,74),
  f_stat = c(81,68,68, 72,81,59),
  author = c(rep("Pulit", 6)),
  doi = c(rep("10.1093/hmg/ddy327", 6)),
  measure = c(rep("weight (kg) / height (m2)", 3), rep("waist circumference (cm) / hip circumference (cm)", 3)),
  adjustment = c(rep("sex (UK Biobank and sex-combined analysis only), age at assessment, age at assessment squared, and assessment centre (UK Biobank only); study specific covariates were included in the GIANT meta-analysis where appropriate",6)),
  transformation = c(rep("inverse rank normal",6)),
  unit = c(rep("normalized SD",6))
)

# cancer table ====
table_cancer <- data.frame(
  data_set = c(rep("cancer", 15)), 
  trait = c("overall","overall","overall",
            "colon", "colon", "colon",
            "distal", "distal", "distal",
            "proximal", "proximal", "proximal",
            "rectal", "rectal", "rectal"),
  sex = c("combined",  "female", "male",
          "combined", "female", "male", 
          "combined", "female", "male", 
          "combined", "female", "male", 
          "combined", "female", "male"),
  
  population = c(rep("European + East Asian", 15)),
  
  N = c("58131",  "6176", "31288",
        "32002", "16105", "15897", 
        "14376", "6689", "7687", 
        "15706", "8520", "7186", 
        "16212", "6541", "9671"),
  N_control = c("67347",  "32820", "34527",
                "64159", "32820", "34527", 
                "64159", "32820", "34527", 
                "64159", "32820", "34527", 
                "64159", "32820", "34527"),
  
  pvalue_threshold = c(rep(5E-8, 15)),
  LD_r2 = c(rep(0.001, 15)),
  LD_window = c(rep(10000, 15)),
  
  SNPs = c(68,26,36,
           44,17,24, 
           31,9,10,
           22,8,3,
           28,13,13),
  f_stat = c(64,55,52,
             57,47,42,
             53,47,48,
             52,47,53,
             63,44,53),
  author = c(rep("Huyghe", 15)),
  doi = c(rep("10.1038/s41588-018-0286-6", 15)),
  measure = c(rep("physician diagnosed", 15)),
  adjustment = c(rep(NA,15)),
  transformation = c(rep(NA,15)),
  unit = c(rep(NA,15))
)


# ferkingstad proteins table ====
a <- fread("data/000_protein_data/exposure_data/000_clumped.txt", stringsAsFactors = T)
a$trait <- paste0(a$sequence_id, "_", a$exposure)
a$trait <- as.factor(a$trait)
a$f_stats <- (a$beta.exposure / a$se.exposure)^2 
f_stats <- a %>%
  group_by(trait) %>%
  summarise(f_stat = mean(f_stats))

n_max <- a %>%
  group_by(trait) %>%
  summarise(N = max(samplesize.exposure))

snp <- a %>%
  group_by(trait) %>%
  summarise(SNPs = length(exposure))

proteins <- left_join(f_stats, n_max)
proteins <- left_join(proteins, snp)
proteins$sex <- "combined"
proteins$N_control <- NA
proteins$data_set <- "proteins"
proteins$population <- "European"
proteins$pvalue_threshold <- 1.8E-9
proteins$LD_r2 <- 0.001
proteins$LD_window <- 10000
proteins$author <- "Ferkingstad"
proteins$doi <- "10.1038/s41588-021-00978-w"
proteins$measure <- "EDTA plasma using SomaScan® by SomaLogic (version 4; 4,907 aptamers)"
proteins$adjustment <- "age, sex, sample age "
proteins$transformation <- "inverse rank normal"
proteins$unit <- "normalized SD"
table_decode <- proteins[,c("data_set","trait","sex",
                            "population","N", "N_control", 
                            "pvalue_threshold", "LD_r2", "LD_window", 
                            "SNPs","f_stat","author","doi",
                            "measure", "adjustment", "transformation", "unit")]

# master table ====
table <- rbind(table_adiposity, table_cancer, table_decode)

write.table(table, "analysis/tables/descriptives.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# protein info ====
proteins <- readxl::read_xlsx("/data/protein_GWAS_ferkingstad_EU_2021/files/protein_info.xlsx")
write.table(proteins, "analysis/tables/proteins.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
