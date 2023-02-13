# combine all adiposity protein outcome files

cd ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data

cat proteins_adiposity/*txt > proteins_adiposity.txt

# add header

head /data/GWAS_data/files/adiposity_GWAS/pulit_EU_2018/processed/bmi_combined.txt

# add header
head -1 /data/GWAS_data/files/adiposity_GWAS/pulit_EU_2018/processed/bmi_combined.txt  > ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/header # get header
cat ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/proteins_adiposity.txt >> ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/header # append data onto header file
mv ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/header ~/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/proteins_adiposity.txt
