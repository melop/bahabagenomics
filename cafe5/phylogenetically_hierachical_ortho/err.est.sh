CAFE=/data/software/CAFE5/bin/cafe5
CPU=24

$CAFE -c $CPU -t tree.tre -i orthogroup.counts.csv -p -e > err.est.log 2>&1 
