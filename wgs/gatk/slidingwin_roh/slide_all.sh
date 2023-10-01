#sSampleList="/data/projects/rcui/bahaha_assembly/wgs/psmc/called_separately/samples.txt"
sSampleList="/data/projects/rcui/bahaha_assembly/wgs/psmc/called_separately/samples_Q1.txt"
arrSamples=( `cut -d" " -f1 $sSampleList` )
arrCovUpper=(  `cut -d" " -f4 $sSampleList` )
arrCovLower=(  `cut -d" " -f3 $sSampleList` )

for nSample in "${!arrSamples[@]}"; do
        sSample=${arrSamples[$nSample]};
        nCovUpper=${arrCovUpper[$nSample]};
        nCovLower=${arrCovLower[$nSample]};
	i=`ls ../$sSample.genotyped*.gz`
	echo sliding_het_finescale.php -i $i -o $sSample.het.win.txt.gz -l $nCovLower -u $nCovUpper
	php sliding_het_finescale.php -i $i -o $sSample.het.win.txt.gz -l $nCovLower -u $nCovUpper  &
done
wait
