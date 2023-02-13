#!/bin/bash
#SBATCH -J obtain_adiposity_instruments
#SBATCH -A mr_adiposity_proteins_crc
#SBATCH --mem=100M
#SBATCH -t 1:00:00

cd ~/001_projects/adiposity_proteins_colorectal_cancer/

## Adiposity measures pulit ====
ls ../000_datasets/adiposity_GWAS/pulit_EU_2018/processed/whr_*.txt | while read f; do awk -F" " 'NR==1{print;next}$9<5e-09' ${f} > ${f}_pulit_snps.txt; done;
ls ../000_datasets/adiposity_GWAS/pulit_EU_2018/processed/bmi_*.txt | while read f; do awk -F" " 'NR==1{print;next}$9<5e-09' ${f} > ${f}_pulit_snps.txt; done;

mv ../000_datasets/adiposity_GWAS/pulit_EU_2018/processed/*_snps.txt data/exposure_data/
