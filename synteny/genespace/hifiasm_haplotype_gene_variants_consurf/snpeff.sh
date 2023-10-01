module load java
module load samtools

sSNPEFF="/public/software/snpEff/exec/snpeff";
ref=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
gff=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/annotations_orthofinder_symbols/btp.longest.addedUPhO.genesymbol2.spgeneid.gff3

mkdir -p snpeff_data/genomes
mkdir -p snpeff_data/ref

ln -sf `realpath $ref` snpeff_data/genomes/ref.fa
ln -sf `realpath $gff` snpeff_data/ref/genes.gff



echo "data.dir = ./snpeff_data/" > snpeff.config
echo "ref.genome : ref" >> snpeff.config



#$sSNPEFF build -c snpeff.config -dataDir `pwd`/snpeff_data/  -gff3 -noCheckProtein -noCheckCds -v ref > ./snpeff.build.log 2>&1
#$sSNPEFF ann -c snpeff.config ref 'merged.coverfilter.vcf.gz' | bgzip -c > 'ann.merged.vcf.gz'

$sSNPEFF ann -c snpeff.config ref 'merged.vcf.gz' | bgzip -c > 'ann.merged.vcf.gz'
tabix ann.merged.vcf.gz

