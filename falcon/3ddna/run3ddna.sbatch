#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 40

Contigs=/data/projects/rcui/bahaha_assembly/falcon/4-polish/cns-output/polished_p_ctgs.fasta
MNF=/data/projects/rcui/bahaha_assembly/falcon/juicer/aligned/merged_nodups.txt
MAPQ=30
GAPSIZE=1000


/data/software/3d-dna-201008/run-asm-pipeline.sh -m haploid -i 100 -r 1 -q $MAPQ --sort-output -g $GAPSIZE $Contigs $MNF
