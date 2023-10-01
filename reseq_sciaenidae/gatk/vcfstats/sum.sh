echo -n "" > summed.stats
for i in *.stats; do
	sB=${i/.genotyped.g.vcf.gz.stats/}
	echo -n "$sB " >> summed.stats
	cat $i | grep -A500 '# DP' | tail -n+3 | sed 's/ +/\t/' | awk 'BEGIN {accp=0;maxdp=0;dp=0;percent25=0;percent975=0} {if (maxdp<$6) {maxdp=$6; dp=$3}; accp+=$7; if (accp<10) {percent25=$3}; if (accp<=97.5) {percent975=$3} } END{print dp" "percent25" "percent975}' >> summed.stats
done
