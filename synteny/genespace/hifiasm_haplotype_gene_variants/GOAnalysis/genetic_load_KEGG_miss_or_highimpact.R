setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/GOanalysis")

library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

datPGID2Ref <- read.table("../pgid2zebrafishhuman.txt", header = T, stringsAsFactors = F)
datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)

datP <- datPresence[,5:ncol(datPresence)]
datH <- datHighImpact[,5:ncol(datHighImpact)]

nSamples <- ncol(datP)
arrSamples <- colnames(datP)
datRet <- datH

for(nIdx in 1:ncol(datP) ) {
  datRet[, nIdx] <- as.integer(datP[ , nIdx] >=0.5 & datH[, nIdx]==0)
}


datOnlyCompHighImpact <- cbind(datPresence[,1:3], datRet);


datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datPGID2Ref, by.x = "pgID", by.y="pgid" , all.x = T, all.y =F)

datZebrafishTrans2GeneID <- as.data.frame( org.Dr.egENSEMBLTRANS)
colnames(datZebrafishTrans2GeneID)[1] <- "ZebrafishGeneID";

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datZebrafishTrans2GeneID, by.x="Zebrafish", by.y = "trans_id", all.x =T, all.y =F);



datHumanTrans2GeneID <- as.data.frame( org.Hs.egENSEMBLTRANS)
colnames(datHumanTrans2GeneID)[1] <- "HumanGeneID";


datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datHumanTrans2GeneID, by.x="Human", by.y = "trans_id", all.x =T, all.y =F);

datHumanGOMap <- as.data.frame(org.Hs.egPATH)
datZebrafishGOMap <- as.data.frame(org.Dr.egPATH)

datMap <- datOnlyCompHighImpact[, c('pgID', "ZebrafishGeneID")]
datMap <- datMap[complete.cases(datMap) , ];

datMap <- merge(datMap, datZebrafishGOMap[,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datMap <- datMap[complete.cases(datMap) , c(3,2)];

datMapHs <- datOnlyCompHighImpact[, c('pgID', "HumanGeneID")]
datMapHs <- datMapHs[complete.cases(datMapHs) , ];

datMapHs <- merge(datMapHs, datHumanGOMap[ ,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datMapHs <- datMapHs[complete.cases(datMapHs) , c(3,2)];

datMap <- rbind(datMap, datMapHs);

datMap <- datMap[!duplicated(datMap),];

colnames(datMap) <- c("term", "gene");


#get gene symbol map
datHumanSymbolMap <- as.data.frame(org.Hs.egSYMBOL)
datZebrafishSymbolMap <- as.data.frame(org.Dr.egSYMBOL)

datSymbol <- datOnlyCompHighImpact[, c('pgID', "ZebrafishGeneID")]
datSymbol <- datSymbol[complete.cases(datSymbol) , ];
datSymbol <- merge(datSymbol, datZebrafishSymbolMap[,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datSymbol <- datSymbol[complete.cases(datSymbol) , c(3,2)];

datSymbolHs <- datOnlyCompHighImpact[, c('pgID', "HumanGeneID")]
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , ];
datSymbolHs <- merge(datSymbolHs, datHumanSymbolMap[,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , c(3,2)];

datSymbol <- rbind(datSymbol , datSymbolHs[ !(datSymbolHs$pgID %in% datSymbol$pgID),] );


fnAddGeneSymbols <- function(datOut) {
  arrSymbolStr <- c();
  if (nrow(datOut) ==0) {
    return(datOut)
  }
  for(nRow in 1:nrow(datOut)) {
    
    arrPGIDs <- unlist(strsplit(datOut$geneID[nRow] ,'/'));
    arrGeneSymbols <- datSymbol[datSymbol$pgID %in% arrPGIDs, 1];
    sStr <- "";
    if (length(arrGeneSymbols)>0) {
      sStr <- paste(arrGeneSymbols, collapse = '/');
    }
    arrSymbolStr <- c(arrSymbolStr, sStr);
  }
  datOut$genesymbols <- arrSymbolStr;
  return(datOut)
}
#get gene symbol

for(sHap in arrSamples) {
  datThis <- datOnlyCompHighImpact[, c('pgID', sHap)];
  arrTargetGenes <- datThis[datThis[,2]!=1 , 'pgID'];
  oEnrichRet <- enricher(
    as.character(arrTargetGenes),
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH",
    as.character(datThis$pgID),
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 1,
    TERM2GENE = datMap ,
    TERM2NAME = NA
  )
  sOut <- paste0("miss_or_highimpact_", sHap,"_KEGG.txt");
  datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=0.5, ]
  write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep="\t", col.names = T, row.names = F);
}

datThis <- datOnlyCompHighImpact[, c('pgID', arrSamples)];
arrTargetGenes <- datThis[rowSums(datThis[,arrSamples])<nSamples , 'pgID'];
oEnrichRet <- enricher(
  as.character(arrTargetGenes),
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  as.character(datThis$pgID),
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = NA
)
sOut <- paste0("miss_or_highimpact_allhaps_KEGG.txt");
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=0.5, ]
write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep="\t", col.names = T, row.names = F);
#viewKEGG(oEnrichRet)
