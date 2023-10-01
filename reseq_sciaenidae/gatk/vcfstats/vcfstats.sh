for i in ../*.genotyped*.gz; do
	sB=`basename $i`
	bcftools stats $i > $sB.stats 2>$sB.log &
done
wait
