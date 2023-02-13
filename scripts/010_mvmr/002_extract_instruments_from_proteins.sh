#!/bin/bash

#SBATCH --job-name=extract-instruments-proteins1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# export
export SNPS=~/001_projects/adiposity_proteins_colorectal_cancer/analysis/008_mvmr/000_unique_snp_list.txt
export PROTEIN_DATA=/data/protein_GWAS_ferkingstad_EU_2021/files/processed
export OUT=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/000_mvmr_data/proteins
export FILES=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/000_mvmr_data/000_file_names_for_mvmr.txt

cd /data/protein_GWAS_ferkingstad_EU_2021/files/processed/

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
for file in `cat $FILES`; do
	zgrep -w -F -f ${SNPS} ${PROTEIN_DATA}/${file} > ${tmp} &&
	mv -- ${tmp} ${OUT}/${file}.instruments
done
