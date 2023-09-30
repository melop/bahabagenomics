#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate merqury

sRd=/public3/group_crf/home/cuirf/bahaha_assembly/hifiasm/female/bahaha.ccs.fastq
k=21

meryl k=$k threads=64  count $sRd output btp.female.meryl

