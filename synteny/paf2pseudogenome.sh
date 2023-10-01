genome1=$1 #/data/projects/rcui/mop/nextpolish/nextpolish_stlfr2/nextpolish/mop.LG.nextdenovo.sgs.2.fa
genome2=$2 #/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.sorted.fa
sp1=$3
sp2=$4

bcftools view -H -i 'TYPE!="snp"' out.$sp1.$sp2.vcf.gz | awk '{ len=length($4); start=$2-1; end=start+len;  print $1"\t"start"\t"end }' > indels.mask.$sp1.$sp2.bed
cat covered.regions.$sp1.$sp2.txt | awk '{print $1"\t"($2-1)"\t"$3}' > covered.regions.$sp1.$sp2.bed
cut -f1-2 $genome1.fai | sort -k1,1 -k2,2n > genome.chrlen.$sp1.txt
bedtools complement -i  covered.regions.$sp1.$sp2.bed -g genome.chrlen.$sp1.txt >  uncovered.regions.$sp1.$sp2.bed
cat indels.mask.$sp1.$sp2.bed uncovered.regions.$sp1.$sp2.bed | sort -k1,1 -k2,2n > masked.regions.$sp1.$sp2.bed
bcftools consensus -i 'TYPE=="snp"'  -f $genome1 out.$sp1.$sp2.vcf.gz > _pseudogenome.$sp1.$sp2.fa
bedtools maskfasta -mc N -fi _pseudogenome.$sp1.$sp2.fa  -fo pseudogenome.$sp1.$sp2.fa -bed  masked.regions.$sp1.$sp2.bed
samtools faidx pseudogenome.$sp1.$sp2.fa
rm _pseudogenome.$sp1.$sp2.fa
