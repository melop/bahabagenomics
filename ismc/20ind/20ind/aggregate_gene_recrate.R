setwd("/data/projects/rcui/bahaha_assembly/ismc/20ind/20ind")
#datRec <- read.table("gene_recrate.txt", header=F)
#datRec <- read.table("gene_recrate.excludeROH.txt", header=F)
datRec <- read.table("gene_recrate.transid.excludeROH.txt", header=F)


datAggRec <- aggregate(datRec[, 8 ], by=list(OrthoID=datRec[,4]), FUN=mean)
colnames(datAggRec)[2] <- "recrate";

datAggRecMinCoord <- aggregate(datRec[, 2 ], by=list(OrthoID=datRec[,4]), FUN=min)
colnames(datAggRecMinCoord)[2] <- "start";

datAggRecMaxCoord <- aggregate(datRec[, 3 ], by=list(OrthoID=datRec[,4]), FUN=max)
colnames(datAggRecMaxCoord)[2] <- "end";

datAggRecScf <- aggregate(datRec[, 1 ], by=list(OrthoID=datRec[,4]), FUN=function(x) {x[1]})
colnames(datAggRecScf)[2] <- "scf";

datAgg <- merge(datAggRecScf, datAggRecMinCoord, by="OrthoID");
datAgg <- merge(datAgg, datAggRecMaxCoord, by="OrthoID");
datAgg <- merge(datAgg, datAggRec, by="OrthoID");

#write.table(datAgg, file="gene_recrate.aggregated.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(datAgg, file="gene_recrate.aggregated.excludeROH.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(datAgg, file="gene_recrate.transid.aggregated.excludeROH.txt", quote = F, sep = "\t", row.names = F, col.names = T)

