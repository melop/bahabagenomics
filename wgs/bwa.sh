CPU=6
REF=/data/projects/rcui/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
for sR1 in ~/../shareddata/bahaha/wgs/trim/BTP*.paired_1.fq.gz; do
        sBase=`basename $sR1`
        sBase=${sBase/.paired_1.fq.gz/}
	sR2=${sR1/.paired_1.fq.gz/.paired_2.fq.gz}
	echo $sBase $sR1 $sR2
        ( bwa mem -t $CPU $REF $sR1 $sR2 \
        | samtools view  -u - | samtools sort - -m 10g -o $sBase.bam && \
	samtools view -h $sBase.bam | awk '{if ($7=="=" && $9<600 && $9>=-600) print $0}'
        ) > $sBase.log 2>&1 &
done
wait
