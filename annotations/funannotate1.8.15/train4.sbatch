#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

source genomedef.sh

CPU=$SLURM_CPUS_PER_TASK

funannotate train -i $GENOME -o trained  \
--left /public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate/hap_M0/alltissues.1.fq.gz \
--right /public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate/hap_M0/alltissues.2.fq.gz  \
--no_trimmomatic \
--no_normalize_reads \
--pasa_min_avg_per_id 85 \
--jaccard_clip \
--memory 100G \
--trinity "$TRINITY_TRANSCRIPT" \
--species "$SPECIES" \
--pasa_db mysql \
--max_intronlen 200000 \
--aligners minimap2 minimap2 \
--strain "$STRAIN" --cpus $CPU 
