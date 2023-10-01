setwd("C:\\data\\BaiduSyncdisk\\projects\\bahaha\\ismc")
datRet <-NULL;
pdf(file="rho.plot.pdf", width=10,height=5)
sink(file = "rho.screen.out.txt")
for (nChr in 10:10) {
  sF <- paste0(nChr, "/out.rho.100kb.bedgraph");
  if (!file.exists(sF)) {
    next;
  }
  dat20ind <- read.table(sF, header=T, sep="\t", fill=T, stringsAsFactors = F)
  dat20ind$chromStart <- dat20ind$chromEnd - 9999;
  dat20ind$chromMid <- (dat20ind$chromStart + dat20ind$chromEnd)/2;
  dat19ind <- dat20ind[, colnames(dat20ind)!="BTP_Q0"] ; #EXCLUDE REF GENOME INDIVIDUAL
  dat19ind$sample_mean <- rowMeans(dat19ind[, 4:(ncol(dat19ind)-3)])
  
  
  
  arrFemales <- grep("BTP_F", colnames(dat19ind));
  arrMales <- grep("BTP_M", colnames(dat19ind));
  
  
  dat19ind$male_sample_mean <- rowMeans(dat19ind[, arrMales])
  dat19ind$female_sample_mean <- rowMeans(dat19ind[, arrFemales])
  
  dat19ind$malefemale_chisqp <-1;
  for(nRow in 1:nrow(dat19ind)) {
    oTest <- t.test(as.numeric(dat19ind[nRow, arrMales]), as.numeric(dat19ind[nRow, arrFemales]) , paired = F);
    dat19ind$malefemale_chisqp[nRow] <- oTest$p.value;
  }
  dat19ind$malefemale_chisqp_fdr <- p.adjust(dat19ind$malefemale_chisqp, method = "fdr")
  
  plot(dat19ind$chromMid, dat19ind$male_sample_mean, type = 'l', col="blue", xlab = "pos", ylab = "rho", main = nChr)
  lines(dat19ind$chromMid, dat19ind$female_sample_mean, col='red')
  
  arrSamples <- colnames(dat19ind)[4:(ncol(dat19ind)-3)];
  
  arrR2 <- c();
  for(sSample in arrSamples) {
    oS <- summary(lm(dat19ind$sample_mean ~ dat19ind[,sSample]))
    #print(oS)
    arrR2 <- c(arrR2, oS$r.squared);
  }
  
  datRet <- rbind(datRet, dat19ind[, c("chrom","chromStart" ,"chromEnd", "sample_mean")]);
  
  cat("#Rsquare", nChr, "\t" ,mean(arrR2), sd(arrR2),"\n");
}

write.table(datRet, file="rho.txt", sep="\t", quote = F, row.names = F, col.names = F)
dev.off();
sink()