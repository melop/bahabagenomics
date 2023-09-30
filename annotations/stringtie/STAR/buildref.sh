#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 32

module load star/2.7.9a

STAR=STAR


source ../genomedef.sh

mkdir -p btpref

$STAR --runThreadN $SLURM_CPUS_PER_TASK \
--genomeSAindexNbases 13 \
--runMode genomeGenerate \
--genomeDir btpref \
--genomeFastaFiles $GENOME

