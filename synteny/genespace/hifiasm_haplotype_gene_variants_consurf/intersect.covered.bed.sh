out=covered.regions.bed

if [ -f $out ]; then rm $out; fi

for i in `ls pairwise.aln.ref.larimichthys_crocea_jimeiU.paf.covered.bed pairwise.aln.ref.Miichthys_miiuyi.paf.covered.bed`; do
	if [ ! -f $out ]; then
		cp $i $out;
	else
		bedtools intersect -a $out -b $i > $out.tmp
		sort -V -k1,1 $out.tmp >  $out
		rm $out.tmp
	fi
done 
