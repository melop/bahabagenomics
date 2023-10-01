setwd("/data/projects/rcui/bahaha_assembly/ismc/20ind/20ind");
options(scipen=999)

#load ROH segments:
arrHomFiles <- Sys.glob("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_roh/*.hom.mid.win.tsv");
#remove M06 since it is same individual as Q1
arrHomFiles <- arrHomFiles[!grepl("BTP_Q1", arrHomFiles)]
length(arrHomFiles)
datROHSegments <- NULL;
for (sF in arrHomFiles) {
  datThis <- read.table(sF, header=T, sep="\t");
  datThis$Group <-  gsub('.hom.mid.win.tsv','', basename(sF), fixed=T);
  datROHSegments <- rbind(datROHSegments , datThis);
}




unlink("rho.excludeROH.txt")

for (nChr in 1:24) {
  sF <- paste0(nChr, "/out.rho.10kb.bedgraph");
  if (!file.exists(sF)) {
    next;
  }
  dat20ind <- read.table(sF, header=T, sep="\t", fill=T, stringsAsFactors = F)
  dat20ind$chromStart <- dat20ind$chromEnd - 9999;
  dat20ind$chromMid <- (dat20ind$chromStart + dat20ind$chromEnd)/2;
  dat19ind <- dat20ind[, colnames(dat20ind)!="BTP_Q0"] ; #EXCLUDE REF GENOME INDIVIDUAL

  dat19ind$chrom <- sChr <- paste0("bahahascf_",nChr);
  
  #set individual ROH segments to NA
  arrInd <- unique(datROHSegments$Group);
  arrInd <- arrInd[arrInd!="BTP_Q0"]
  for(sInd in arrInd  ) {
    datROHThis <- datROHSegments[datROHSegments$Group==sInd & datROHSegments$chr==sChr,]
    if (nrow(datROHThis) == 0) {
      next;
    }
    for(nROHRow in 1:nrow(datROHThis)) {
      dat19ind[dat19ind$chromMid >= datROHThis[nROHRow, 'homstart'] & dat19ind$chromMid <= datROHThis[nROHRow, 'homend'],  sInd ] <-NA;
    }
  }
  
  
  
  dat19ind$sample_mean <- rowMeans(dat19ind[,arrInd], na.rm = T)

  write.table(dat19ind[, c("chrom","chromStart" ,"chromEnd", "sample_mean")], file="rho.excludeROH.txt", sep="\t", quote = F, row.names = F, col.names = F,append = T)
  

}


