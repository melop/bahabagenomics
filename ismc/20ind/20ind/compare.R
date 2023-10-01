
options(scipen=999)

unlink("rho.txt")
pdf(file="rho.plot.pdf", width=10,height=5)
sink(file = "rho.screen.out.txt")
for (nChr in 1:24) {
  sF <- paste0(nChr, "/out.rho.10kb.bedgraph");
  if (!file.exists(sF)) {
    next;
  }
  dat20ind <- read.table(sF, header=T, sep="\t", fill=T, stringsAsFactors = F)
  dat20ind$chromStart <- dat20ind$chromEnd - 9999;
  dat20ind$chromMid <- (dat20ind$chromStart + dat20ind$chromEnd)/2;
  dat19ind <- dat20ind[, colnames(dat20ind)!="BTP_Q0"] ; #EXCLUDE REF GENOME INDIVIDUAL
  dat19ind$sample_mean <- rowMeans(dat19ind[, 4:(ncol(dat19ind)-3)])
  
  dat19ind$chrom <- paste0("bahahascf_",nChr);
  
  plot(dat19ind$chromMid, dat19ind$sample_mean, type = 'l', col="blue", xlab = "pos", ylab = "rho", main = nChr)
  
  arrSamples <- colnames(dat19ind)[4:(ncol(dat19ind)-3)];
  
  arrR2 <- c();
  for(sSample in arrSamples) {
    oS <- summary(lm(dat19ind$sample_mean ~ dat19ind[,sSample]))
    #print(oS)
    arrR2 <- c(arrR2, oS$r.squared);
  }
  
  write.table(dat19ind[, c("chrom","chromStart" ,"chromEnd", "sample_mean")], file="rho.txt", sep="\t", quote = F, row.names = F, col.names = F,append = T)
  
  cat("#Rsquare", nChr, "\t" ,mean(arrR2), sd(arrR2),"\n");
}

dev.off();
sink()
