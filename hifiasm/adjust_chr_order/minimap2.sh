#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 64

module load minimap2
hap0=../3ddna/female.hap1.p_ctg.0.review.fasta
hap1=../3ddna_hap2/female.hap2.p_ctg.0.review.fasta

falconhap0=/public3/group_crf/home/cuirf/bahaha_assembly/annotations/repeatmask/female.hap0/adjust_chr_order/hap.F0.fa
falconhap1=/public3/group_crf/home/cuirf/bahaha_assembly/annotations/repeatmask/female.hap1/adjust_chr_order/hap.F1.fa


minimap2 -t 64  -x asm5 $hap0  $falconhap0  > aln.0.to.hap0.paf
minimap2 -t 64  -x asm5 $hap0  $falconhap1  > aln.0.to.hap1.paf

minimap2 -t 64  -x asm5 $hap1  $falconhap0  > aln.1.to.hap0.paf
minimap2 -t 64  -x asm5 $hap1  $falconhap1  > aln.1.to.hap1.paf
