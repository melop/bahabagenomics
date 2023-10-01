sSNPEFF="/public/software/snpEff/exec/snpeff";

#no need to filter out a priori, php script performs filter
#outgroup_cov_bed=/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf/covered.regions.bed

#bcftools view -R $outgroup_cov_bed -O z merged.vcf.gz  >  merged.coverfilter.vcf.gz
#tabix merged.coverfilter.vcf.gz


$sSNPEFF ann -c snpeff.config ref 'merged.vcf.gz' | bgzip -c > 'merged.snpeff.vcf.gz'
tabix merged.snpeff.vcf.gz



