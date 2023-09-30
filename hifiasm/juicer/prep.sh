#!/usr/bin/bash
#SBATCH -p blade,gpu,himem
#SBATCH -c 1

module load bwa
module load python3
module load samtools

genome=../female.hap1.p_ctg.fa

genome=`realpath $genome`

mkdir -p references
ln -sf $genome references/ref.fa
bwa index references/ref.fa
samtools faidx references/ref.fa

mkdir -p restriction_sites
cut -f1,2 references/ref.fa.fai > restriction_sites/ref.chrom.sizes

cd restriction_sites
python /public3/group_crf/software/juicer-1.6/misc/generate_site_positions.py DpnII ref ../references/ref.fa

cd ..
mkdir fastq
cd fastq
ln -s /public3/group_crf/home/cuirf/bahaha_assembly/hifiasm/female/E100033076_L01_33.paired_1.fq.gz ./read_R1.fastq.gz
ln -s /public3/group_crf/home/cuirf/bahaha_assembly/hifiasm/female/E100033076_L01_33.paired_2.fq.gz ./read_R2.fastq.gz
