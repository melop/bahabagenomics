#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 1

module load java

Contigs=/public3/group_crf/home/cuirf/bahaha_assembly/hifiasm/female/female.hap1.p_ctg.fa
MNF=/public3/group_crf/home/cuirf/bahaha_assembly/hifiasm/female/juicer/aligned/merged_nodups.txt
MAPQ=30
GAPSIZE=1000


/public/software/3d-dna-201008/run-asm-pipeline.sh -m haploid -i 100 -r 1 -q $MAPQ --sort-output -g $GAPSIZE $Contigs $MNF
