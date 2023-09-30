#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 36
#SBATCH --mem=30G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate merqury

hap0=/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate/hifiasm_hap_F0/genome.rename.fa
hap1=/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate/hifiasm_hap_F1/genome.rename.fa

ln -sf $hap0 ./hap0.fa
ln -sf $hap1 ./hap1.fa
merqury.sh btp.female.meryl hap0.fa hap1.fa out



