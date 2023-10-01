#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 24


module load samtools
module load clusterbasics
module load bcftools
module load java
module load minimap2

echo run snpeff_short_alleles.php
php snpeff_short_alleles.php
