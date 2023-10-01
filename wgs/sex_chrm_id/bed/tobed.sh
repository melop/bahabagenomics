cut -f1,2 /data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa.fai | awk '{print $1"\t"($2+100001)}' > chrom.sizes

echo 'track type=bedGraph name="difference between female and male depth" description="difference between female and male depth" priority=20' > m.f.dp_diff.bed
cat ../m.f.dp_diff.plot.tsv | tail -n+2 | awk '{ print $1"\t"($2-1)"\t"$3"\t"$(NF);}' | sort -k1,1 -k2,2n | awk '{mean=($2+$3)/2; printf "%s\t%d\t%d\t%f\n", $1, (mean-1), mean, $4}' >> m.f.dp_diff.bed
bedGraphToBigWig m.f.dp_diff.bed chrom.sizes m.f.dp_diff.bw
#bgzip m.f.dp_diff.bed
#tabix m.f.dp_diff.bed.gz


echo 'track type=bedGraph name="difference between female and male heterozygosity" description="difference between female and male heterozygosity" priority=20' > m.f.het_diff.bed
cat ../m.f.het_diff.plot.tsv | tail -n+2 | awk '{ print $1"\t"($2-1)"\t"$3"\t"$(NF);}' | sort -k1,1 -k2,2n | awk '{mean=($2+$3)/2; printf "%s\t%d\t%d\t%f\n", $1, (mean-1), mean, $4}' >> m.f.het_diff.bed
bedGraphToBigWig m.f.het_diff.bed chrom.sizes m.f.het_diff.bw
#bgzip m.f.het_diff.bed
#tabix m.f.het_diff.bed.gz

echo 'track type=bedGraph name="difference between female and male depth chr10" description="difference between female and male depth chr10" priority=20' > m.f.dp_diff.5kb.chr10.bed
cat ../m.f.dp_diff.5kb.chr10.plot.tsv | tail -n+2 | awk '{ print $1"\t"($2-1)"\t"$3"\t"$(NF);}' | sort -k1,1 -k2,2n | awk '{mean=($2+$3)/2; printf "%s\t%d\t%d\t%f\n", $1, (mean-1), mean, $4}' >> m.f.dp_diff.5kb.chr10.bed
bedGraphToBigWig m.f.dp_diff.5kb.chr10.bed chrom.sizes m.f.dp_diff.5kb.chr10.bw
#bgzip m.f.dp_diff.5kb.chr10.bed
#tabix m.f.dp_diff.5kb.chr10.bed.gz

echo 'track type=bedGraph name="FST female-male" description="FST female-male" priority=20' > m.f.fst.bed
cat ../fst_female_vs_male.windowed.weir.fst | tail -n+2 | awk '{ print $1"\t"($2-1)"\t"$3"\t"$5;}' | sort -k1,1 -k2,2n | awk '{mean=($2+$3)/2; printf "%s\t%d\t%d\t%f\n", $1, (mean-1), mean, $4}' >> m.f.fst.bed
bedGraphToBigWig  m.f.fst.bed chrom.sizes  m.f.fst.bw
#bgzip m.f.fst.bed
#tabix m.f.fst.bed.gz

cat ../markers_XY.vcf | awk '{print $1"\t"($2-1)"\t"$2}' | sort -k1,1 -k2,2n  > markers_XY.bed
bgzip markers_XY.bed
tabix markers_XY.bed.gz

cat ../markers_ZW.vcf | awk '{print $1"\t"($2-1)"\t"$2}' | sort -k1,1 -k2,2n  > markers_ZW.bed
bgzip markers_ZW.bed
tabix markers_ZW.bed.gz
