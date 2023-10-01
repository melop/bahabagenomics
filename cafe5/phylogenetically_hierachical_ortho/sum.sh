dir=runresults_k4
dir=runresults_k5
dir=runresults
dir=runresults_k2
outF=$dir.sum.txt
echo -n "" > $outF
for i in `seq 1 20`; do
	lnl=`grep 'Final -lnL:' $dir/$i/est.log | tail -n1 | cut -d' ' -f3`
	echo "$i	$lnl" >> $outF
done

cat $outF | sort -k2,2n > _tmp.txt
mv _tmp.txt $outF
