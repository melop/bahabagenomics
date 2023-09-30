#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 32


source genomedef.sh

module load samtools
module load minimap2
module load stringtie

mkfifo all.fulllength.ont.cdna.fa
samtools  fasta  ~/../shareddata/bahaha/rna/PB/RT-1018/bc1018/split.bc1018_5p--bc1018_3p.subreads.bam >  all.fulllength.ont.cdna.fa &
REF=$GENOME
CPU=$SLURM_CPUS_PER_TASK

minimap2 -ax splice -t $CPU --cs -u b -G 1000000 $REF all.fulllength.ont.cdna.fa | samtools sort --reference $REF -@ 12 -o btp.ont.cdna.alignments.bam - > aln.log 2>&1

stringtie  --mix  -p $CPU -o longreads_assm.out.gtf STAR/alltissues.bam btp.ont.cdna.alignments.bam > stringtie.log 2>&1
/public/software/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl longreads_assm.out.gtf $REF > btp.assembled_transcripts.fa
