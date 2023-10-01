bcftools view -O z -s ^BTP_Q1 ../joinedcalled.snps.indels.vcf.gz > exclude.btpq1.vcf.gz
bcftools index exclude.btpq1.vcf.gz
bcftools merge --missing-to-ref  -o merged.vcf.gz -O z exclude.btpq1.vcf.gz  out*vcf.gz
rm exclude.btpq1.vcf.gz
