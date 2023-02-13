#!/bin/bash

#SBATCH --job-name=protein-instruments-filelist14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /data/protein_GWAS_ferkingstad_EU_2021/files/processed/

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
for file in `cat filelist14`; do
    gzip -d -c ${file} > ${file}.unzipped
    awk -F" " 'NR==1{print;next}$8<5e-08' ${file}.unzipped > ${tmp} &&
    mv -- ${tmp} ${file}.snps
    rm ${file}.unzipped
    mv ${file}.snps /data/protein_GWAS_ferkingstad_EU_2021/files/instruments/
done

