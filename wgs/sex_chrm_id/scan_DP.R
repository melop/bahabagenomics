setwd("~/bahaha_assembly/wgs/sex_chrm_id/");


datDP <- NULL;
arrFiles<- Sys.glob("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_DP_100kb_10kb/*.DP.win.txt")

for( sF in arrFiles) {
  sFileName <- basename(sF);
  datThis <- read.table(sF, header = F, sep="\t", stringsAsFactors = T, na.strings = 'NAN')
  sSample <- gsub(".DP.win.txt", "", sFileName);
  datThisDP <- data.frame(scf=datThis$V1, winstart=datThis$V2,winend=datThis$V3, winmid=datThis$V4, DP=datThis$V5)
  colnames(datThisDP)[5] <- sSample;
  if (is.null(datDP)) {
    datDP <- datThisDP;
  } else {
    datDP <- merge(datDP, datThisDP[,c(1,2,5)], by = c('scf', 'winstart') )
  }
  
}

datDP$chrnum <-as.integer( gsub("[^0-9]", "", datDP$scf))
datDP <- datDP[ order(datDP$chrnum, datDP$winstart), ];
datDP$chrnum <-NULL

arrFemales <- grep("_F" , colnames(datDP))
arrMales <- grep("_M" , colnames(datDP))


datDP$ttest.p <- -1.0;
datDP$female_avg <- 0.0;
datDP$male_avg <- 0.0;
datDP$dpdiff <- 0.0;
scf <- "";
arrFemalesMeanDP <- colMeans(datDP[, arrFemales] , na.rm = T)
arrMalesMeanDP <- colMeans(datDP[, arrMales] , na.rm = T)

for(i in 1:nrow(datDP)) {
  arrFVals <- datDP[i, arrFemales] / arrFemalesMeanDP;
  arrMVals <- datDP[i, arrMales] / arrMalesMeanDP;
  arrFVals <- arrFVals[!is.na(arrFVals)]
  arrMVals <- arrMVals[!is.na(arrMVals)]
  sScfThis <- as.character(datDP[i, 'scf']);
  if (sScfThis!=scf) {
    print(sScfThis)
    scf <- sScfThis
  }
  
  if (length(arrMVals) < 3 || length(arrFVals)<3 ) {
    datDP[i, 'ttest.p'] <- NA;
    datDP[i, 'dpdiff'] <- NA;
  } else {
    oT <- t.test(arrFVals, arrMVals, paired=F);
    datDP[i, 'ttest.p'] <- oT$p.value;
    datDP[i, 'female_avg'] <- oT$estimate[1];
    datDP[i, 'male_avg'] <- oT$estimate[2];
    datDP[i, 'dpdiff'] <- datDP[i, 'female_avg'] - datDP[i, 'male_avg'] 
  }
}

hist(log10(datDP$ttest.p) );

datDPNoNA <- datDP[!is.na(datDP$ttest.p) , ]


datDPNoNA$log10_ttest.p <- -log10(datDPNoNA$ttest.p)
datDPNoNA[datDPNoNA$dpdiff < 0 , 'log10_ttest.p'] <- log10(datDPNoNA[datDPNoNA$dpdiff < 0 , 'ttest.p'])

write.table(datDPNoNA, file = "m.f.dp_diff.plot.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
pdf(file="m.f.dp_diff.plot.pdf", width = 20, height = 5)
for(nChr in 1:24) {
  datMiss10 <- datDPNoNA[datDPNoNA$scf==paste0('bahahascf_', nChr),]
  plot(datMiss10$winmid, datMiss10$log10_ttest.p, type = 'l', lwd=1, col="red", ylim = c(-20,20), xlab=paste("chromosome", nChr), ylab="log10(p)")

}
dev.off();


