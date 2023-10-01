sga=/data/software/sga/src/SGA/sga
sp=$1
r1=$2
r2=$3
nSubSampleLn=500000000 #about 45x coverage
CPU=$4




mkdir -p $sp

(
mkfifo $sp/r1.fastq
mkfifo $sp/r2.fastq

zcat -f $r1 | head -n $nSubSampleLn > $sp/r1.fastq &
zcat -f $r2 | head -n $nSubSampleLn > $sp/r2.fastq &


$sga preprocess --pe-mode 1  $sp/r1.fastq $sp/r2.fastq > $sp/$sp.fastq
cd $sp
$sga index -a ropebwt --no-reverse -t $CPU  $sp.fastq
$sga preqc -v  -t $CPU  $sp.fastq >  $sp.preqc && ( rm *.fastq ; rm *.bwt ; rm *.sai )
#sga-preqc-report.py mygenome.preqc sga/src/examples/*.preqc

) >> $sp/log.txt 2>&1
