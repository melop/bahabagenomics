cat repeats.bed | cut -f4 | cut -d'/' -f1 | sort | uniq -c > repeat.type.count.txt
