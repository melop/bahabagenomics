grep HIGH SV.w.indelmodifier.polarized.out.txt  | grep -P "2\t2\t2" | grep -vP "\t1" | grep -vP "NA"  > fixed.highimpact.txt
