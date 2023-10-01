setwd("~/bahaha_assembly/wgs/sex_chrm_id/");

datMissing <- NULL;
datHet <- NULL;
arrFiles<- Sys.glob("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_het_100kb_10kb/*.het.win.txt")

for( sF in arrFiles) {
  sFileName <- basename(sF);
  datThis <- read.table(sF, header = F, sep="\t", stringsAsFactors = T, na.strings = 'NAN')
  sSample <- gsub(".het.win.txt", "", sFileName);
  datThisMiss <- data.frame(scf=datThis$V1, winstart=datThis$V2,winend=datThis$V3, winmid=datThis$V4, missing=datThis$V8)
  colnames(datThisMiss)[5] <- sSample;
  if (is.null(datMissing)) {
    datMissing <- datThisMiss;
  } else {
    datMissing <- merge(datMissing, datThisMiss[,c(1,2,5)], by = c('scf', 'winstart') )
  }
  
  datThisHet <- data.frame(scf=datThis$V1, winstart=datThis$V2,winend=datThis$V3, winmid=datThis$V4, het=datThis$V9)
  colnames(datThisHet)[5] <- sSample;
  if (is.null(datHet)) {
    datHet <- datThisHet;
  } else {
    datHet <- merge(datHet, datThisHet[,c(1,2,5)], by = c('scf', 'winstart') )
  }
  
}

datMissing$chrnum <-as.integer( gsub("[^0-9]", "", datMissing$scf))
datMissing <- datMissing[ order(datMissing$chrnum, datMissing$winstart), ];
datMissing$chrnum <-NULL

datHet$chrnum <-as.integer( gsub("[^0-9]", "", datHet$scf))
datHet <- datHet[ order(datHet$chrnum, datHet$winstart), ];
datHet$chrnum <-NULL

arrFemales <- grep("_F" , colnames(datHet))
arrMales <- grep("_M" , colnames(datHet))


datHet$ttest.p <- -1.0;
datHet$female_avg <- 0.0;
datHet$male_avg <- 0.0;
datHet$hetdiff <- 0.0;
scf <- "";
for(i in 1:nrow(datHet)) {
  arrFVals <- datHet[i, arrFemales];
  arrMVals <- datHet[i, arrMales];
  arrFVals <- arrFVals[!is.na(arrFVals)]
  arrMVals <- arrMVals[!is.na(arrMVals)]
  sScfThis <- as.character(datHet[i, 'scf']);
  if (sScfThis!=scf) {
    print(sScfThis)
    scf <- sScfThis
  }
  
  if (length(arrMVals) < 3 || length(arrFVals)<3 ) {
    datHet[i, 'ttest.p'] <- NA;
    datHet[i, 'hetdiff'] <- NA;
  } else {
    oT <- t.test(arrFVals, arrMVals, paired=F);
    datHet[i, 'ttest.p'] <- oT$p.value;
    datHet[i, 'female_avg'] <- oT$estimate[1];
    datHet[i, 'male_avg'] <- oT$estimate[2];
    datHet[i, 'hetdiff'] <- datHet[i, 'female_avg'] - datHet[i, 'male_avg'] 
  }
}

hist(log10(datHet$ttest.p) );

datHetNoNA <- datHet[!is.na(datHet$ttest.p) , ]
#View(datHetNoNA[datHetNoNA$ttest.p<1e-5,])

#quantile(datHetNoNA$ttest.p, c(0.005, 0.995))

datHetNoNA$log10_ttest.p <- -log10(datHetNoNA$ttest.p)
datHetNoNA[datHetNoNA$hetdiff < 0 , 'log10_ttest.p'] <- log10(datHetNoNA[datHetNoNA$hetdiff < 0 , 'ttest.p'])

write.table(datHetNoNA, file = "m.f.het_diff.plot.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
pdf(file="m.f.het_diff.plot.pdf", width = 20, height = 5)
for(nChr in 1:24) {
datHet10 <- datHetNoNA[datHetNoNA$scf==paste0('bahahascf_', nChr),]
plot(datHet10$winmid, datHet10$log10_ttest.p, type = 'l', lwd=1, col="blue", ylim = c(-10,10), xlab=paste("chromosome", nChr), ylab="log10(p)")
}
dev.off();






arrFemales <- grep("_F" , colnames(datMissing))
arrMales <- grep("_M" , colnames(datMissing))


datMissing$ttest.p <- -1.0;
datMissing$female_avg <- 0.0;
datMissing$male_avg <- 0.0;
datMissing$missdiff <- 0.0;
scf <- "";
arrFemalesMeanMiss <- colMeans(datMissing[, arrFemales] , na.rm = T)
arrMalesMeanMiss <- colMeans(datMissing[, arrMales] , na.rm = T)

for(i in 1:nrow(datMissing)) {
  arrFVals <- datMissing[i, arrFemales] / arrFemalesMeanMiss;
  arrMVals <- datMissing[i, arrMales] / arrMalesMeanMiss;
  arrFVals <- arrFVals[!is.na(arrFVals)]
  arrMVals <- arrMVals[!is.na(arrMVals)]
  sScfThis <- as.character(datMissing[i, 'scf']);
  if (sScfThis!=scf) {
    print(sScfThis)
    scf <- sScfThis
  }
  
  if (length(arrMVals) < 3 || length(arrFVals)<3 ) {
    datMissing[i, 'ttest.p'] <- NA;
    datMissing[i, 'missdiff'] <- NA;
  } else {
    oT <- t.test(arrFVals, arrMVals, paired=F);
    datMissing[i, 'ttest.p'] <- oT$p.value;
    datMissing[i, 'female_avg'] <- oT$estimate[1];
    datMissing[i, 'male_avg'] <- oT$estimate[2];
    datMissing[i, 'missdiff'] <- datMissing[i, 'female_avg'] - datMissing[i, 'male_avg'] 
  }
}

hist(log10(datMissing$ttest.p) );

datMissNoNA <- datMissing[!is.na(datMissing$ttest.p) , ]


datMissNoNA$log10_ttest.p <- -log10(datMissNoNA$ttest.p)
datMissNoNA[datMissNoNA$missdiff < 0 , 'log10_ttest.p'] <- log10(datMissNoNA[datMissNoNA$missdiff < 0 , 'ttest.p'])

write.table(datMissNoNA, file = "m.f.miss_diff.plot.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
pdf(file="m.f.het.miss_diff.plot.pdf", width = 20, height = 5)
for(nChr in 1:24) {
  datMiss10 <- datMissNoNA[datMissNoNA$scf==paste0('bahahascf_', nChr),]
  plot(datMiss10$winmid, datMiss10$log10_ttest.p, type = 'l', lwd=1, col="red", ylim = c(-20,20), xlab=paste("chromosome", nChr), ylab="log10(p)")
  datHet10 <- datHetNoNA[datHetNoNA$scf==paste0('bahahascf_', nChr),]
  lines(datHet10$winmid, datHet10$log10_ttest.p, type = 'l', lwd=1, col="blue", ylim = c(-20,20) )
  
}
dev.off();


