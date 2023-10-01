#bedtools intersect -a /data/projects/rcui/bahaha_assembly/hyphy/exportAln_genespace/groupid2coord.txt -b rho.txt -wa -wb > gene_recrate.txt
#bedtools intersect -a /data/projects/rcui/bahaha_assembly/hyphy/exportAln_genespace/groupid2coord.txt -b rho.excludeROH.txt -wa -wb > gene_recrate.excludeROH.txt

bedtools intersect -a /data/projects/rcui/bahaha_assembly/hyphy/exportAln_genespace/transid2coord.txt -b rho.excludeROH.txt -wa -wb > gene_recrate.transid.excludeROH.txt
