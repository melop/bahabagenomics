#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=50G

source genomedef.sh

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

CPU=$SLURM_CPUS_PER_TASK

funannotate update -i trained --cpus $CPU --pasa_db mysql --alt_transcripts 0.3
