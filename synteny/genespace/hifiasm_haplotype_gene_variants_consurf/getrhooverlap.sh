cat consurf.polarized.out.txt | awk '{ start=$2-1;  print $1"\t"start"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}'  > consurf_af.bed
bedtools intersect -a consurf_af.bed -b ~/bahaha_assembly/ismc/20ind/20ind/rho.excludeROH.txt -wa -wb   > consurf_af.rec.bed
