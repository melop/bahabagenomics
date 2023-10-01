
CAFE=/data/software/CAFE5/bin/cafe5
CPU=32

#for i in `seq 1 10`; do
for i in `seq 11 20`; do
	mkdir -p runresults/$i
$CAFE -c $CPU -t tree.tre -i orthogroup.counts.tsv -p -eresults/Base_error_model.txt -k 3 -o  runresults/$i >  runresults/$i/est.log 2>&1 
done
