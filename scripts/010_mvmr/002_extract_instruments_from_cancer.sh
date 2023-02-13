#!/bin/bash

#SBATCH --job-name=extract-instruments-cancer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# export
export SNPS=~/001_projects/adiposity_proteins_colorectal_cancer/analysis/008_mvmr/000_unique_snp_list.txt
export CANCER_DATA=/data/GWAS_data/files/huyghe_2018_PMID30510241/processed
export OUT=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/000_mvmr_data/cancer

cd /data/GWAS_data/files/huyghe_2018_PMID30510241/processed/
ls *annotated.* > filenames.txt
export FILES=/data/GWAS_data/files/huyghe_2018_PMID30510241/processed/filenames.txt

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
for file in `cat $FILES`; do
	grep -f ${SNPS} ${CANCER_DATA}/${file} > ${tmp} &&
	mv -- ${tmp} ${OUT}/${file}.instruments
done


