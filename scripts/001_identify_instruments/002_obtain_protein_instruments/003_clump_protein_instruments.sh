#!/bin/bash

#SBATCH --job-name=protein-instruments-filelist1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /home/leem/001_projects/adiposity_proteins_colorectal_cancer

Rscript scripts/001_identify_instruments/002_obtain_protein_instruments/003_clump_protein_instruments.R
