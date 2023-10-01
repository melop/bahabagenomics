#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 24


module load diamond
module load clusterbasics
module load R

php get_length_outliers_isoforms.php
