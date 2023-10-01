for i in *_pseudogenome2.fasta; do
 cat $i | grep -io 'N' | wc -l > $i.countN.txt &
done
wait
