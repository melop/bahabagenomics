gvcf=/data/projects/rcui/bahaha_assembly/wgs/gatk/all.samples.genotyped.vcf.gz
maxDP=889
maxMissing=0.5
winSize=100000
winStep=10000
bcftools view -O z -i 'TYPE=="snp" && FILTER!="LowQual"' $gvcf > var.only.vcf.gz
tabix var.only.vcf.gz
vcftools --gzvcf var.only.vcf.gz --weir-fst-pop females.txt --weir-fst-pop males.txt --out fst_female_vs_male --min-alleles 2  --max-alleles 2 --fst-window-size $winSize --fst-window-step $winStep --max-meanDP $maxDP --max-missing $maxMissing
