module load minimap2
module load bcftools
module load samtools

genome1=$1 #/data/projects/rcui/mop/nextpolish/nextpolish_stlfr2/nextpolish/mop.LG.nextdenovo.sgs.2.fa
genome2=$2 #/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.sorted.fa
sp1=$3
sp2=$4
algorithm=$5

paftools=paftools.js

minimap2 $algorithm  -t6   --secondary=no   --cs $genome1 $genome2 | sort -k6,6 -k8,8n --parallel=4 > pairwise.aln.$sp1.$sp2.paf

cut -f1,3,4,5,6,8,9,12 pairwise.aln.$sp1.$sp2.paf | awk '$8>=60 {print}' > pairwise.simp.$sp1.$sp2.txt

( $paftools call -s $sp2 -q 60 -f $genome1 -L10000 -l1000 pairwise.aln.$sp1.$sp2.paf | bgzip -c | bcftools norm -w 99999999 -O z -f $genome1  - > out.$sp1.$sp2.vcf.gz ) 2>out.vcf.info
tabix out.$sp1.$sp2.vcf.gz

