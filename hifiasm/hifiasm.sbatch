#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 64

module load hifiasm


hifiasm -o bahaba_female -t 64 --h1 E100033076_L01_33.paired_1.fq.gz --h2 E100033076_L01_33.paired_2.fq.gz bahaha.ccs.fastq
