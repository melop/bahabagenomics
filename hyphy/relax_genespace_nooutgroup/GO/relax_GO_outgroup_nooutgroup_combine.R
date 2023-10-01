setwd("/data/projects/rcui/bahaha_assembly/hyphy/relax_genespace_nooutgroup/GO");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

nFDRCutoff <- 0.01
nRawPCutoff <- 1
bRelaxed <- T; #True  - relaxed genes, False - intensified genes
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

datRelax <- read.table("../sum_BTP.txt", header=F, sep="\t", fill = T, quote = "")
datRelax2 <- read.table("../../relax_genespace/sum_BTP.txt", header=F, sep="\t", fill = T, quote = "")

datRelaxM <- merge(datRelax[, c(1,5,9)], datRelax2[, c(1,5,9)], by="V1", all.x=T, all.y=T)

datRelaxM$V5 <- 1;
datRelaxM$V9 <- 0;

for(i in 1:nrow(datRelaxM)) {
  if (is.na(datRelaxM[i , 'V5.x'])) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.y'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.y'];
  } else if (is.na(datRelaxM[i , 'V5.y'])) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.x'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.x'];
  } else if (datRelaxM[i , 'V5.y'] < datRelaxM[i , 'V5.x'] ) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.y'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.y'];   
  } else {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.x'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.x'];    
  }
}

datRelax <- datRelaxM[, c(1, 6, 7)];
datRelax$fdr <- p.adjust(datRelaxM$V5, method = "fdr")
datGroupID2TransID <- read.table("../../exportAln_genespace/groupid2transid.BTP.txt", sep="\t", header=T)
colnames(datGroupID2TransID) <- c("OrthoID", 'Gene');
datRelax <- merge(datRelax, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)

write.table(datRelax, file="relax.outnooutcombined.txt", sep="\t", col.names=T, row.names=F, quote=F);

quit(); #pause

datOnlyCompHighImpact <- datRNAID2ZebraFishHuman; #merge(datConsurf, datRNAID2ZebraFishHuman, by.x = "Gene", by.y=1 , all.x = T, all.y =F)


datZebrafishTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/zebrafish_ncbi_ensembl_id_map.txt", sep="\t", header=F) ;#as.data.frame( org.Dr.egENSEMBLTRANS)
datZebrafishTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datZebrafishTrans2GeneID$V5)
datZebrafishTrans2GeneID <- datZebrafishTrans2GeneID[, c(2,5)]
colnames(datZebrafishTrans2GeneID) <- c("ZebrafishGeneID", 'trans_id');

# datZebrafishTrans2GeneID <- as.data.frame( org.Dr.egENSEMBLTRANS)
# colnames(datZebrafishTrans2GeneID)[1] <- "ZebrafishGeneID";

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datZebrafishTrans2GeneID, by.x="Zebrafish", by.y = "trans_id", all.x =T, all.y =F);


datHumanTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/human_ncbi_ensembl_id_map.txt", sep="\t", header=F); #as.data.frame( org.Hs.egENSEMBLTRANS)
datHumanTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datHumanTrans2GeneID$V5)
datHumanTrans2GeneID <- datHumanTrans2GeneID[, c(2,5)]
colnames(datHumanTrans2GeneID) <- c("HumanGeneID", 'trans_id');

# datHumanTrans2GeneID <- as.data.frame( org.Hs.egENSEMBLTRANS)
# colnames(datHumanTrans2GeneID)[1] <- "HumanGeneID";

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

arrIsRelaxed <- datRelax$V9 < 1;
sRelaxed <- "relaxed";
if ( ! bRelaxed ) {
  arrIsRelaxed <- datRelax$V9 > 1;
  sRelaxed <- "intensified";
}
arrTargetGenes <- unique(datRelax[ datRelax$V5<= nRawPCutoff & datRelax$fdr <= nFDRCutoff & arrIsRelaxed , 'Gene']);
oEnrichRet <- enricher(
  as.character(arrTargetGenes),
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  as.character(unique(datRelax[ , 'Gene'])),
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)

sOut <- paste0("outAndNoOutCombinded.", sRelaxed,"_fdr_", nFDRCutoff, "_p_", nRawPCutoff, "_GO_",sCategory,".txt");
datOut <- oEnrichRet@result[ oEnrichRet@result$pvalue <= 0.05 , ]
View(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep="\t", col.names = T, row.names = F);
#}


###volcano plot
library(ggplot2)
library(ggrepel)
library(dplyr)
datRelaxPlot <- datRelax;

# nMinFdr <- min(datRelaxPlot$fdr[datRelaxPlot$fdr>0])
# datRelaxPlot$fdr[datRelaxPlot$fdr>0] <- nMinFdr;

datRelaxPlot$P <- datRelaxPlot$V5
mMinP <- min(datRelaxPlot$P[datRelaxPlot$P>0]);
datRelaxPlot$P[datRelaxPlot$P < mMinP] <- mMinP;

nMaxK <- max(datRelaxPlot$V9);
datRelaxPlot$V9[datRelaxPlot$V9 < (1/nMaxK)] <-  (1/nMaxK);

datRelaxPlot$logK <- log10(datRelaxPlot$V9);

nFDR001P <- max(datRelaxPlot[datRelaxPlot$fdr <=0.001, 'P'])
nFDR01P <- max(datRelaxPlot[datRelaxPlot$fdr <=0.01, 'P'])
nFDR05P <- max(datRelaxPlot[datRelaxPlot$fdr <=0.05, 'P'])

datRelaxPlot <- datRelaxPlot %>% 
  mutate(
    Significance = case_when(
      P <= nFDR05P & P > nFDR01P ~ "fdr 0.05", 
      P <= nFDR01P & P > nFDR001P ~ "fdr 0.01",
      P <= nFDR001P ~ "fdr 0.001", 
      TRUE ~ "Unchanged")
  )


p3 <- ggplot(datRelaxPlot, aes(logK, -log10(P))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[10]*"K")) + 
  ylab(expression("-log"[10]*"P")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 


datRelaxPlot$Selection <- "Relaxed"
datRelaxPlot$Selection[datRelaxPlot$logK >0] <- "Intensified"

top <- 20
top_genes <- bind_rows(
  datRelaxPlot %>% 
    filter(Selection == 'Intensified' & fdr<=0.05) %>% 
    arrange(fdr, desc(abs(logK))) %>% 
    head(top),
  datRelaxPlot %>% 
    filter(Selection == 'Relaxed' & fdr<=0.05) %>% 
    arrange(fdr, desc(abs(logK))) %>% 
    head(top)
)

colnames(top_genes)[5] <- "geneID"
top_genes <- fnAddGeneSymbols(top_genes)
top_genes$genesymbols[top_genes$genesymbols==""] <- ""; #top_genes$geneID[top_genes$genesymbols==""]
top_genes$genesymbols <- str_to_lower(top_genes$genesymbols);

p3 <- p3 + geom_label_repel(data = top_genes, max.overlaps = 100, 
                   mapping = aes(logK, -log10(P), label = genesymbols),
                   size = 4)

datPosSel <- read.table("../../bustedmh_genespace/sumBTP.txt", header=F, sep="\t", quote = "", fill = T);
arrPosGenes <- datPosSel$V1[datPosSel$V9<0.05]
datPosSel <- read.table("../../bustedmh_genespace_nooutgroup/sumBTP.txt", header=F, sep="\t", quote = "", fill = T);
arrPosGenes <- unique(c(arrPosGenes, datPosSel$V1[datPosSel$V9<0.05]));

top_pos_genes <- datRelaxPlot[datRelaxPlot$V1 %in% arrPosGenes, ]
colnames(top_pos_genes)[5] <- "geneID"
top_pos_genes <- fnAddGeneSymbols(top_pos_genes)
top_pos_genes$genesymbols[top_pos_genes$genesymbols==""] <- ""; #top_genes$geneID[top_genes$genesymbols==""]
top_pos_genes$genesymbols <- str_to_lower(top_pos_genes$genesymbols);


p3 <- p3 + geom_label_repel( data = top_pos_genes, max.overlaps = 100, fontface = 'bold',
                            mapping = aes(logK, -log10(P), label = genesymbols),
                            size = 4)


p3

ggsave("volcanoplot.pdf", width=6,height=5)

nIntenseGenesSig <- nrow(datRelax[datRelax$fdr<0.05 & datRelax$V9>1, ])
nRelaxGenesSig <- nrow(datRelax[datRelax$fdr<0.05 & datRelax$V9<1, ])
nIntenseGenesNonSig <- nrow(datRelax[datRelax$fdr>0.05 & datRelax$V9>1, ])
nRelaxGenesNonSig <- nrow(datRelax[datRelax$fdr>0.05 & datRelax$V9<1, ])

matTest <- matrix(c(nIntenseGenesSig, nIntenseGenesNonSig, nRelaxGenesSig, nRelaxGenesNonSig), nrow = 2 )
fisher.test(matTest)
