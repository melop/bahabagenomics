cat ~/bahaha_assembly/release1.0/annotations_orthofinder_symbols/btp.longest.addedUPhO.genesymbol2.spgeneid.gff3  | awk '{if ($3=="CDS") print $1"\t"($4-1)"\t"$5}'  > cds.bed
bedtools subtract -A -a cds.bed -b repeats.bed > cds.norepeats.bed
