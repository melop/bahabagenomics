#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 1

JUICERDIR=/data/software/juicer

$JUICERDIR/scripts/juicer.sh \
-D $JUICERDIR \
-g ref \
-d `pwd` \
-z `pwd`/references/ref.fa \
-p `pwd`/restriction_sites/ref.chrom.sizes \
-y `pwd`/restriction_sites/ref_DpnII.txt \
-q long -l long \
> screen.log 2>&1

