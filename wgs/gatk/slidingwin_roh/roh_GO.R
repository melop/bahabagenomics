setwd("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_roh");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

bRelaxed <- F; #True  - relaxed genes, False - intensified genes
sRefSp <- "BahahaTaipingensis";
sCategory <- 'BP' ; #biological process
datPGID2Ref <- read.table("/data/projects/rcui/bahaha_assembly/synteny/genespace/bahaha_haplotypes_zebrafish_human_forGO/synorthos.txt", header = T, stringsAsFactors = F)

datRNAID2ZebraFishHuman <- datPGID2Ref[, c(sRefSp, paste0(sRefSp,'.1'), 'Human', 'Human.1', 'Zebrafish', 'Zebrafish.1') ];
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, sRefSp]) , sRefSp ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, sRefSp]) , paste0(sRefSp,'.1') ]
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Human']) , 'Human' ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, "Human"]) , 'Human.1' ]
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) , 'Zebrafish' ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) , 'Zebrafish.1' ]
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[, c(sRefSp, 'Human', 'Zebrafish') ]
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!is.na(datRNAID2ZebraFishHuman[,sRefSp]) , ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!(is.na(datRNAID2ZebraFishHuman[, 'Human']) & is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) ), ];
datRNAID2ZebraFishHuman$Human <- sub(";.*", "", datRNAID2ZebraFishHuman$Human)
datRNAID2ZebraFishHuman$Zebrafish <- sub(";.*", "", datRNAID2ZebraFishHuman$Zebrafish)

colnames(datRNAID2ZebraFishHuman)[1] <- "Gene";
datRNAID2ZebraFishHuman$Gene <- gsub('.', '_', datRNAID2ZebraFishHuman$Gene, fixed = T);


datOnlyCompHighImpact <- datRNAID2ZebraFishHuman; #merge(datConsurf, datRNAID2ZebraFishHuman, by.x = "Gene", by.y=1 , all.x = T, all.y =F)


datZebrafishTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/zebrafish_ncbi_ensembl_id_map.txt", sep="\t", header=F) ;#as.data.frame( org.Dr.egENSEMBLTRANS)
datZebrafishTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datZebrafishTrans2GeneID$V5)
datZebrafishTrans2GeneID <- datZebrafishTrans2GeneID[, c(2,5)]
colnames(datZebrafishTrans2GeneID) <- c("ZebrafishGeneID", 'trans_id');

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datZebrafishTrans2GeneID, by.x="Zebrafish", by.y = "trans_id", all.x =T, all.y =F);


datHumanTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/human_ncbi_ensembl_id_map.txt", sep="\t", header=F); #as.data.frame( org.Hs.egENSEMBLTRANS)
datHumanTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datHumanTrans2GeneID$V5)
datHumanTrans2GeneID <- datHumanTrans2GeneID[, c(2,5)]
colnames(datHumanTrans2GeneID) <- c("HumanGeneID", 'trans_id');


datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datHumanTrans2GeneID, by.x="Human", by.y = "trans_id", all.x =T, all.y =F);

datHumanGOMap <- as.data.frame(org.Hs.egGO)
datZebrafishGOMap <- as.data.frame(org.Dr.egGO)


datMap <- datOnlyCompHighImpact[, c('Gene', "ZebrafishGeneID")]
datMap <- datMap[complete.cases(datMap) , ];

datMap <- merge(datMap, datZebrafishGOMap[datZebrafishGOMap$Ontology == sCategory,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datMap <- datMap[complete.cases(datMap) , c(3,2)];

datMapHs <- datOnlyCompHighImpact[, c('Gene', "HumanGeneID")]
datMapHs <- datMapHs[complete.cases(datMapHs) , ];

datMapHs <- merge(datMapHs, datHumanGOMap[datHumanGOMap$Ontology ==sCategory ,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datMapHs <- datMapHs[complete.cases(datMapHs) , c(3,2)];

datMap <- rbind(datMap, datMapHs);

datMap <- datMap[!duplicated(datMap),];

colnames(datMap) <- c("term", "gene");

#get gene symbol map
datHumanSymbolMap <- as.data.frame(org.Hs.egSYMBOL)
datZebrafishSymbolMap <- as.data.frame(org.Dr.egSYMBOL)

datSymbol <- datOnlyCompHighImpact[, c('Gene', "ZebrafishGeneID")]
datSymbol <- datSymbol[complete.cases(datSymbol) , ];
datSymbol <- merge(datSymbol, datZebrafishSymbolMap[,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datSymbol <- datSymbol[complete.cases(datSymbol) , c(3,2)];

datSymbolHs <- datOnlyCompHighImpact[, c('Gene', "HumanGeneID")]
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , ];
datSymbolHs <- merge(datSymbolHs, datHumanSymbolMap[,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , c(3,2)];

datSymbol <- rbind(datSymbol , datSymbolHs[ !(datSymbolHs$Gene %in% datSymbol$Gene),] );
datSymbol <- datSymbol[!duplicated(datSymbol),];

fnAddGeneSymbols <- function(datOut) {
  arrSymbolStr <- c();
  if (nrow(datOut) ==0) {
    return(datOut)
  }
  for(nRow in 1:nrow(datOut)) {
    
    arrPGIDs <- unlist(strsplit(datOut$geneID[nRow] ,'/'));
    arrGeneSymbols <- datSymbol[datSymbol$Gene %in% arrPGIDs, 1];
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

#for(sHap in arrSamples) {

datGeneCoord <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/mRNA.coord.txt", sep="\t", header=T, quote="", fill = T)

arrMidpoints <- (datGeneCoord$start + datGeneCoord$end)/2



arrHomFiles <- Sys.glob("*.hom.mid.win.tsv");
#remove Q0 because this individual is much more inbred.
arrHomFiles <- arrHomFiles[!grepl("BTP_Q", arrHomFiles)]
length(arrHomFiles)
datROHSegments <- NULL;
arrColors <- as.character(wes_palette("Zissou1", length(arrHomFiles), type = "continuous"))
nColor <- 1;
for (sF in arrHomFiles) {
  datThis <- read.table(sF, header=T, sep="\t");
  datThis$Group <- gsub('BTP_','', gsub('.hom.mid.win.tsv','', sF, fixed=T), fixed=T);
  datThis$Colors <- arrColors[nColor]
  nColor <- nColor + 1;
  datROHSegments <- rbind(datROHSegments , datThis);
}

arrIsInRange <- rep(F, nrow(datGeneCoord))
for (sChr in unique(datROHSegments$chr) ) {
  datROHSegmentsChr <- datROHSegments[datROHSegments$chr==sChr,];
  for(nRow in 1:nrow(datROHSegmentsChr)) {
    arrIsInRange <- (arrIsInRange | (datGeneCoord$chr==sChr & arrMidpoints>=datROHSegmentsChr[nRow,'homstart'] & arrMidpoints<=datROHSegmentsChr[nRow,'homend']));
  }
}


arrTargetGenes <- unique(datGeneCoord[ arrIsInRange , 'mRNAId']);

oEnrichRet <- enricher(
  as.character(arrTargetGenes),
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  as.character(unique(datGeneCoord[ , 'mRNAId'])),
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)

        
sOut <- paste0("rho.GO.tsv.txt");
datOut <- oEnrichRet@result[ oEnrichRet@result$pvalue <=0.05, ]
#View(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep="\t", col.names = T, row.names = F);
#}


