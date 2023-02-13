#!/bin/bash

#SBATCH --job-name=extract-MVMR-instruments
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:0:00
#SBATCH --mem=100000M

# export
export SNPS=~/001_projects/adiposity_proteins_colorectal_cancer/analysis/008_mvmr/000_unique_snp_list.txt

export bmi_combined=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/bmi_combined.txt
export bmi_males=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/bmi_male.txt
export bmi_females=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/bmi_female.txt

export whr_combined=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/whr_combined.txt
export whr_males=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/whr_male.txt
export whr_females=/data/GWAS_data/files/pulit_2018_PMID30239722/processed/whr_female.txt

export OUT=~/001_projects/adiposity_proteins_colorectal_cancer/data/000_mvmr_data/adiposity

# grep instruments from each GWAS files for BMI
grep -w -F -f ${SNPS} ${bmi_combined} > ${OUT}/bmi_combined.txt
grep -w -F -f ${SNPS} ${bmi_males} > ${OUT}/bmi_males.txt
grep -w -F -f ${SNPS} ${bmi_females} > ${OUT}/bmi_females.txt

grep -w -F -f ${SNPS} ${whr_combined} > ${OUT}/whr_combined.txt
grep -w -F -f ${SNPS} ${whr_males} > ${OUT}/whr_males.txt
grep -w -F -f ${SNPS} ${whr_females} > ${OUT}/whr_females.txt

