#!/usr/bin/bash
#SBATCH -p short
#SBATCH -c 32
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate busco5.3.2

busco -c $SLURM_CPUS_PER_TASK -i ../*.longest_isoform.prot.fa  -m prot -o longest_prot -l actinopterygii_odb10 \
--offline --download_path /public2/shareddatabase/busco/download/ > busco.full.log 2>&1
