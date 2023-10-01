setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf");
datConsurfRec <- read.table("consurf_af.rec.bed", header=F, sep="\t")
datConsurfRec$chr <- as.integer(gsub("bahahascf_", "", datConsurfRec$V1))

datConsurfRec <- datConsurfRec[ rowSums(datConsurfRec[,6:9])<4, ];
datConsurfRec <- datConsurfRec[complete.cases(datConsurfRec[,6:9]),];#only keep segregating variants

datConsurfRec1 <- datConsurfRec[datConsurfRec$V5==1, ];
datConsurfRec9 <- datConsurfRec[datConsurfRec$V5>=9, ];

for(nHap in 6:9) {
  
  arrMutatedRecRate <- datConsurfRec9[datConsurfRec9[,nHap]==1 , 13];
  arrNonMutatedRecRate <- datConsurfRec9[datConsurfRec9[,nHap]==0, 13];
  # print(median(arrMutatedRecRate))
  # print(median(arrNonMutatedRecRate))
  # print(wilcox.test(arrMutatedRecRate, arrNonMutatedRecRate))
  # 
  print(mean(arrMutatedRecRate))
  print(mean(arrNonMutatedRecRate))
  print(t.test(arrMutatedRecRate, arrNonMutatedRecRate))
  
}