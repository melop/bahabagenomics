CPU=64
REF=/data/projects/rcui/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa
READS=/data/projects/rcui/bahaha_assembly/falcon/bahaha.ccs.fastq
sBase=BTP_Q0

( minimap2 -a -x map-pb -t $CPU $REF $READS \
        | samtools view  -u - | samtools sort - -m 10g -o $sBase.bam ) > map.log 2>&1
