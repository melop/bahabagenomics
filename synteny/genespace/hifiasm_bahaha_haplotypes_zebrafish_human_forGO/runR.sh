#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G

source /public/apps/miniconda3/etc/profile.d/conda.sh

conda activate orthofinder
module load mcscanx
module load R/4.2.0

Rscript run.R 
