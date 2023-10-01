#for sIn in *.het.win.txt.gz; do
for sIn in *Q1.het.win.txt.gz; do

        sSample=${sIn/.het.win.txt.gz/};
	echo $sSample
	sOut=${sSample}.hom.mid.win.tsv
	echo Rscript plot_het_mid.R $sIn $sOut
	Rscript plot_het_mid.R $sIn $sOut &
done
wait
