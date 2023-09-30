#!/usr/bin/bash
#SBATCH -p short
#SBATCH -c 48

hap=both
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate busco5.3.2

busco -c $SLURM_CPUS_PER_TASK -i ../hap.reassigned.$hap.fa  -m genome -o hap$hap -l actinopterygii_odb10 \
--offline --download_path /public2/shareddatabase/busco/download/ > busco.hap$hap.log 2>&1
