#!/usr/bin/bash
#SBATCH -p blade,gpu
#SBATCH -c 20
ref=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
for sHap in `ls ~/bahaha_assembly/annotations/protein_refs/larimichthys_crocea_jimeiU/scf.softmasked.fa ~/bahaha_assembly/annotations/other_refs/Miichthys_miiuyi/GCA_001593715.1_ASM159371v1_genomic.fna`; do
	sStem=`dirname $sHap`
	sStem=`basename $sStem`
	bash aln_out.sh $ref $sHap ref $sStem '-x asm20' > $sStem.log 2>&1 &
done
wait
