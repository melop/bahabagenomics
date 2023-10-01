winsize=1000000

cut -f1,2 ~/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa.fai | awk '{if ($2>='$winsize') print $0}'  > chr.def.txt

bedtools makewindows -g chr.def.txt -w $winsize > windows.bed
