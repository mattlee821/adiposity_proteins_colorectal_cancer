#!/bin/bash

#SBATCH --job-name=MR-colorectal-proteins
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd ~/001_projects/adiposity_proteins_colorectal_cancer

Rscript scripts/008_colorectal_proteins/001_colorectal_proteins_MR.R