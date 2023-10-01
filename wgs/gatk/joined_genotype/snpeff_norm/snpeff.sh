#snpeffjar=~/../../software/snpEff/snpEff.jar
REF=Bahaha_taipingensis_1.0
genome=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
bcftools view -O z -M 2 -m 2 -i 'TYPE=="snp" || TYPE=="indel"' ../joincalled.genotyped.g.vcf.gz | \
bcftools norm -w 99999999 -O z -f $genome -  > joinedcalled.snps.indels.vcf.gz
#java -Xmx26g -jar $snpeffjar $REF joinedcalled.snps.indels.vcf.gz | bgzip -c > joinedcalled.snps.indels.snpeff.vcf.gz


