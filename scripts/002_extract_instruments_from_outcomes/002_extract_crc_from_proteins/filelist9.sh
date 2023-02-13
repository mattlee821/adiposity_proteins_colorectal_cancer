#!/bin/bash

#SBATCH --job-name=extract-crc-instruments-proteins_filelist9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# adiposity
export SNP_LIST=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/exposure_data/snp_list_crc.txt
export PROTEIN_DATA=/data/protein_GWAS_ferkingstad_EU_2021/files/processed
export OUTCOMES=/home/leem/001_projects/adiposity_proteins_colorectal_cancer/data/outcome_data/colorectal_proteins/

cd /data/protein_GWAS_ferkingstad_EU_2021/files/processed/

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
for file in `cat filelist9`; do
	zgrep -w -F -f ${SNP_LIST} ${file} > ${tmp} &&
	mv -- ${tmp} ${file}.colorectal_instruments
    mv ${file}.colorectal_instruments ${OUTCOMES}
done
