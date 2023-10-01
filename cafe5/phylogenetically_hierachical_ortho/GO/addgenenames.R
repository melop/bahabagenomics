setwd("/data/projects/rcui/bahaha_assembly/cafe5/phylogenetically_hierachical_ortho/GO")
sF <- "BahahaTaipingensis_expandedfam.GO.tsv";
sF <- "BahahaTaipingensis_shrinkedfam.GO.tsv";

sp <- "BahahaTaipingensis";
datGO <- read.table(sF, header = T, sep="\t")
datGeneSymbolMap <- read.table("/data/projects/rcui/bahaha_assembly/annotations/UPhO_morespp/assigngenesymbols_orthofinder/orthofinder_genesymbols.tsv", header=T, sep="\t", fill=T)
datGO$PossibleGeneSymbols <- NA
datGO$PossibleGeneNames <- NA

# datGeneSymbolMap$AllDescriptions <- URLdecode(datGeneSymbolMap$AllDescriptions)
# datGeneSymbolMap$AllGeneSymbols <- URLdecode(datGeneSymbolMap$AllGeneSymbols)

for (i in 1:nrow(datGO)) {
  arrGeneIDs <- str_trim(unlist(strsplit(x = datGO[i, 'geneID'], split = '/')));
  arrPossibleSymbols <- c();
  arrPossibleNames <- c();
  for(sGeneID in arrGeneIDs) {
    arrWhich <- grep( pattern =sGeneID, datGeneSymbolMap[,sp]); 
    arrPossibleSymbols <- c(arrPossibleSymbols, datGeneSymbolMap[arrWhich, 'AllGeneSymbols']);
    arrPossibleNames <- c(arrPossibleNames, datGeneSymbolMap[arrWhich, 'AllDescriptions']);
  }
  
  arrPossibleSymbols <- unique(str_trim(unlist(strsplit(x = arrPossibleSymbols, split = '\\|'))));
  arrPossibleNames <- unique(str_trim(unlist(strsplit(x = arrPossibleNames, split = '\\|'))));
  datGO$PossibleGeneSymbols[i] <- paste(arrPossibleSymbols, collapse = '|') 
  datGO$PossibleGeneNames[i] <- paste(arrPossibleNames, collapse = ' | ') 
}

datGO$PossibleGeneNames <- URLdecode(datGO$PossibleGeneNames)

write.table(datGO, file=paste0(sF, ".genenames.tsv"), sep = "\t", col.names = T, row.names = F, quote = F );
