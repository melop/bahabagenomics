sp1=$1
sp2=$2

paftools=/data/software/minimap2-2.20/misc/paftools.js

$paftools call -q 60 -L10000 -l1000 pairwise.aln.$sp1.$sp2.paf | grep -P "R\t" | cut -f2-4 > covered.regions.$sp1.$sp2.txt 
