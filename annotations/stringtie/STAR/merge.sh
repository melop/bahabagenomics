#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 1


module load samtools

mkdir -p lnret
for i in  mapped/*/*.bam; do
	s=`dirname $i`
	s=`basename $s`
	ln -sf `realpath $i` lnret/$s.bam

done

samtools merge -r alltissues.bam lnret/*.bam
