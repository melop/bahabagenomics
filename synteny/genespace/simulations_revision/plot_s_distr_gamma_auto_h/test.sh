tail -n+2 simdef.txt | \
while read Chr	TotalGen	SeedingPopSize	PopSize	h	h_sd	s	s_sd	RecPerChr	TopPercentileForBreeding	MeanHomDel	SDHomDel	MeanDel	SDDel	GenotypedMarkerPerc	Reps; do
	echo $Chr $TotalGen
done

