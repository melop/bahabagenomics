CPU=12
REF=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
for i in ../*.bam; do
	sID=`basename $i`
	sID=${sID/.bam/}
	echo $sID
	bash gatk.sh $sID $i $REF $CPU > $sID.main.log 2>&1 &
done
wait
