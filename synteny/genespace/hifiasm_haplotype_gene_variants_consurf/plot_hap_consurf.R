setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf");
datConsurf <- read.table("consurf.polarized.out.txt", sep="\t")

hist(datConsurf$V6)

datHiImpact <- datConsurf[datConsurf$V6>=9, ];
arrNAs <- apply(datHiImpact[, 7:10], 2, FUN= function(x) {sum(as.integer(is.na(x)))} )
arrAdjLoad <- colSums(datHiImpact[, 7:10], na.rm = T) / (1-arrNAs / nrow(datHiImpact))

arrAdjLoad
