#!/bin/bash

#SBATCH --job-name=extract-protein-instruments-adiposity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# adiposity
export SNP_LIST=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/000_protein_data/exposure_data/snp_list.txt
export ADIPOSITY_DATA=/home/leem/001_projects/000_datasets/adiposity_GWAS/pulit_EU_2018/processed/
export OUTCOMES=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/proteins_adiposity

# grep instruments from each CRC GWAS file for adiposity instruments
grep -f ${SNP_LIST} ${ADIPOSITY_DATA}/bmi_male.txt > ${OUTCOMES}/bmi_male.txt
