#!/usr/bin/bash
#SBATCH -p blade,gpu
#SBATCH -c 40
ref=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
for sHap in /public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate/hifiasm_hap_*/genome.rename.fa; do
	sStem=`dirname $sHap`
	sStem=`basename $sStem`
	sStem=${sStem/hifiasm_/}
	sStem=${sStem/_/.}
	echo bash aln.sh $ref $sHap ref $sStem
	bash aln.sh $ref $sHap ref $sStem '-x asm5' > $sStem.log 2>&1 &
done

wait
