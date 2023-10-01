bcftools view -R covered.regions.bed -O z  merged.vcf.gz >  merged.coverfilter.vcf.gz
tabix merged.coverfilter.vcf.gz
