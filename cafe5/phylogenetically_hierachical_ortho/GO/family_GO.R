setwd("/data/projects/rcui/bahaha_assembly/cafe5/phylogenetically_hierachical_ortho/GO");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)


sInDir <- "../runresults/18/";
sRefSp <- sSp <- "BahahaTaipingensis";
nFDRCutoff <- 0.2
sCategory <- 'BP' ; #biological process
datPGID2Ref <- read.table("/data/projects/rcui/bahaha_assembly/synteny/genespace/bahaha_haplotypes_zebrafish_human_forGO/rundir/orthofinder/Results_Aug13/Orthogroups/Orthogroups.tsv", header = T, stringsAsFactors = F, fill=T, sep="\t")
datFam2Gene <- read.table("/data/projects/rcui//bahaha_assembly/synteny/genespace/morespp/rundir/orthofinder/Results_Mar18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene$HOG <- gsub('N0.', '', datFam2Gene$HOG , fixed = T)

datChange <- read.table(paste0(sInDir, "/Gamma_change.tab" ), sep="\t", header=T);
datChange <- datChange[, c(1, grep(sSp, colnames(datChange), fixed = T)) ];

datBranchProb <- read.table(paste0(sInDir, "/Gamma_branch_probabilities.tab" ), sep="\t", header=T, comment.char = '/', fill = T);

datChange <- merge(datChange , datBranchProb[, c(1,grep(sSp, colnames(datBranchProb)))], by.x = "FamilyID", by.y=1)

datChange$fdr <- p.adjust(datChange[, 3])

arrExpandedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] > 0  ,1]
arrShrinkedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] < 0  ,1]
arrControlFam <- datChange[ ,1]

fnFam2Gene <- function(arr) {
  arrGenes <- datFam2Gene[datFam2Gene$HOG %in% arr, sSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

arrExpandedGenes <- fnFam2Gene(arrExpandedFam);
arrShrinkedGenes <- fnFam2Gene(arrShrinkedFam);
arrControlGenes <- fnFam2Gene(arrControlFam);

#prepare mapping tables

datRNAID2ZebraFishHuman <- datPGID2Ref[, c(sRefSp,  'Human',  'Zebrafish') ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[(!is.na(datRNAID2ZebraFishHuman[,sRefSp])) & (!datRNAID2ZebraFishHuman[,sRefSp]=='') , ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!(datRNAID2ZebraFishHuman[, 'Human']=='' & datRNAID2ZebraFishHuman[, 'Zebrafish']=='' ), ];


datZebrafishTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/zebrafish_ncbi_ensembl_id_map.txt", sep="\t", header=F) ;#as.data.frame( org.Dr.egENSEMBLTRANS)
datZebrafishTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datZebrafishTrans2GeneID$V5)
datZebrafishTrans2GeneID <- datZebrafishTrans2GeneID[, c(2,5)]
colnames(datZebrafishTrans2GeneID) <- c("ZebrafishGeneID", 'trans_id');
datHumanTrans2GeneID <- read.table("/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/GO/human_ncbi_ensembl_id_map.txt", sep="\t", header=F); #as.data.frame( org.Hs.egENSEMBLTRANS)
datHumanTrans2GeneID$V5 <- gsub('\\.[0-9]*' , "", datHumanTrans2GeneID$V5)
datHumanTrans2GeneID <- datHumanTrans2GeneID[, c(2,5)]
colnames(datHumanTrans2GeneID) <- c("HumanGeneID", 'trans_id');


datHumanGOMap <- as.data.frame(org.Hs.egGO)
datZebrafishGOMap <- as.data.frame(org.Dr.egGO)

datMap <- NULL;
sTERM2GO <- paste0(sSp,"_term2go.tsv");
bAppend <- F;

sink(file=paste0(sSp, ".gomap.log"))
for(i in 1:nrow(datRNAID2ZebraFishHuman)) {
  arrZebrafishTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 'Zebrafish'], split = ',')));
  arrHumanTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 'Human'], split = ',')));
  arrSpTransID <- str_trim(unlist(strsplit(x = datRNAID2ZebraFishHuman[i, 1], split = ',')));
  arrZebrafishGeneID <- datZebrafishTrans2GeneID[datZebrafishTrans2GeneID$trans_id %in% arrZebrafishTransID, 'ZebrafishGeneID']
  arrHumanGeneID <- datHumanTrans2GeneID[datHumanTrans2GeneID$trans_id %in% arrHumanTransID, 'HumanGeneID']
  
  arrGOs <- datZebrafishGOMap[datZebrafishGOMap$gene_id %in% arrZebrafishGeneID & datZebrafishGOMap$Ontology ==sCategory, 'go_id']
  arrGOs <- c(arrGOs, datHumanGOMap[datHumanGOMap$gene_id %in% arrHumanGeneID & datHumanGOMap$Ontology ==sCategory, 'go_id']);
  arrGOs <- unique(arrGOs)
  if (length(arrGOs) == 0) {
    cat("GO not found: ", datRNAID2ZebraFishHuman[i, 1], "\n");
    next;
  } else {
    cat("GO found: ", datRNAID2ZebraFishHuman[i, 1], "\n");
  }
  for(sSpTransID in arrSpTransID) {
    datToBind <- data.frame(term=arrGOs, gene=sSpTransID, stringsAsFactors = F);
    write.table(datToBind, file = sTERM2GO, append = bAppend, sep = "\t", col.names = (!bAppend), row.names = F, quote = F)
    bAppend<-T;
    
  }
}
sink()

datMap <- read.table(sTERM2GO, sep="\t", header=T, stringsAsFactors = T);


#enrichment test:

# arrExpandedGenes 
# arrShrinkedGenes 
# arrControlGenes

sum(as.integer(unique(datMap$gene) %in% arrExpandedGenes))
sum(as.integer(unique(datMap$gene) %in% arrControlGenes))

#View( go2term(unique(datMap$term[datMap$gene %in% arrExpandedGenes])) )
arrExpandedGenes[!(arrExpandedGenes %in% unique(datMap$gene))]
                                                
oEnrichRet <- enricher(
  gene=arrExpandedGenes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 0.1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
#View(datOut)
write.table(datOut, paste0(sSp, "_expandedfam.GO.tsv"), sep = "\t", col.names = T, row.names = F, quote = F);


oEnrichRet <- enricher(
  gene=arrShrinkedGenes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 0.1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
#View(datOut)
write.table(datOut, paste0(sSp, "_shrinkedfam.GO.tsv"), sep = "\t", col.names = T, row.names = F, quote = F);

arrUnMappedExpanded <- arrExpandedGenes[!(arrExpandedGenes %in% unique(datMap$gene))];

arrShow <-c();
for(sUnmapped in arrUnMappedExpanded) {
  arrShow <- c(arrShow, grep(sUnmapped, datPGID2Ref$BahahaTaipingensis));
}

#View(datPGID2Ref[unique(arrShow),]);
