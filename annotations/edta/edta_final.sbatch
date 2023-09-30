#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=120G
CPU=32
genome=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
repbaselib=/public2/shareddatabase/repbase/RepeatMasker.lib
mu=7.752e-10 # estimated from divergence between BTP and Miichthys, this is per generation, per year is 7.752e-10, assuming 10 year generation time


source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

perl /public/software/conda_envs/EDTA/bin/EDTA.pl  --genome  $genome --step all  --curatedlib $repbaselib --step final --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --u "$mu"  --threads $CPU 
