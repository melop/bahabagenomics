cat SV.w.indelmodifier.polarized.out.txt | awk -F"\t" -v OFS='\t' '{ start=$2-1;  printf("%s",$1"\t"start"\t"$2"\t");$1="";$2="";print}'  > consurf_af_sv.w.indelmodifier.bed
bedtools intersect -a consurf_af_sv.w.indelmodifier.bed -b /public3/group_crf/home/cuirf/bahaha_assembly/ismc/20ind/20ind/rho.excludeROH.txt -wa -wb   > consurf_af_sv.w.indelmodifier.rec.bed
