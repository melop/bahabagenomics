#!/usr/bin/bash
#SBATCH -p long,blade,gpu,himem,hugemem
#SBATCH -c 20


module load clusterbasics

php confirm_missing_genes.php
