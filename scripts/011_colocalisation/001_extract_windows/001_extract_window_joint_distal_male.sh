#!/bin/bash

#SBATCH --job-name=extract-window-jdm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-10:0:00
#SBATCH --mem=100000M

cd ~/001_projects/adiposity_proteins_colorectal_cancer

Rscript scripts/011_colocalisation/001_extract_windows/001_extract_window_joint_distal_male.R
