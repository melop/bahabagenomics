if [ -f roh.tsv ]; then
rm roh.tsv

fi


for sIn in *.hom.mid.win.tsv; do
        sSample=${sIn/.hom.mid.win.tsv/};
	echo $sSample
	Rscript roh_stats.R $sIn $sSample  >> roh.tsv
done
wait
