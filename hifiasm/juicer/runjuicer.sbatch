#!/usr/bin/bash
#SBATCH -p blade,gpu
#SBATCH -c 1
#SBATCH --mem=24G

JUICERDIR=/public3/group_crf/software/juicer-1.6/

$JUICERDIR/scripts/juicer.sh \
-D /public3/group_crf/software/juicer-1.6/ \
-g ref \
-d `pwd` \
-s DpnII \
-z `pwd`/references/ref.fa \
-p `pwd`/restriction_sites/ref.chrom.sizes \
-y `pwd`/restriction_sites/ref_DpnII.txt \
-q blade,gpu,himem \
-l blade,gpu,himem \
-t 8 \
> juicer.log 2>&1
