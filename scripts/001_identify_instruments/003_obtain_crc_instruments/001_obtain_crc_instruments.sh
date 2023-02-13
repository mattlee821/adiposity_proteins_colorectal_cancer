#!/bin/bash
#SBATCH -J obtain_crc_instruments
#SBATCH -A mr_adiposity_proteins_crc
#SBATCH --mem=100M
#SBATCH -t 0-1:00:00

cd ~/001_projects/adiposity_proteins_colorectal_cancer/

## CRC ====
ls /data/GWAS_data/files/huyghe_2018_PMID30510241/processed/*TBL.annotated.txt | while read f; do awk -F" " 'NR==1{print;next}$10<5e-08' ${f} > ${f}snps; done;
ls /data/GWAS_data/files/huyghe_2018_PMID30510241/processed/*tsv.annotated.txt | while read f; do awk -F" " 'NR==1{print;next}$12<5e-08' ${f} > ${f}snps; done;
mv /data/GWAS_data/files/huyghe_2018_PMID30510241/processed/*snps data/exposure_data/
