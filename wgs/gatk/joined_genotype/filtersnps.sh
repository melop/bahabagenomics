
bcftools view -O z  -v "snps"  -e 'SOR>2 || QD<5 || MQRankSum<-12.5 || INFO/DP>1000 || INFO/DP<228  || ExcessHet>30' joincalled.genotyped.g.vcf.gz > joincalled.genotyped.g.filtered.snps.vcf.gz
