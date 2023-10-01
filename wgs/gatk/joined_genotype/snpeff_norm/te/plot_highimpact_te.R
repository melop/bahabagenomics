setwd("~/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/te")

pdf(file="te_correlate.AF.pdf", width = 5, height=5)

datHighImpactVsTE <- read.table("cov.perwin.bed", header=F, sep="\t")
colnames(datHighImpactVsTE) <- c('chr', 'start', 'end', 'te_fraction' , 'dna_te_fraction', 'ltr_fraction', 'cds_fraction', 'af0_20', 'af20_40', 'af40_60', 'af60_80', 'af80_99', 'af99_100')

# nMinCDSFraction <- quantile(datHighImpactVsTE$cds_fraction, c(0.025))
# datHighImpactVsTE <- datHighImpactVsTE[datHighImpactVsTE$cds_fraction >=nMinCDSFraction,];

bExcludeZeros <- T;

#plot(datHighImpactVsTE[datHighImpactVsTE$chr=="bahahascf_9", 'start'],datHighImpactVsTE[datHighImpactVsTE$chr=="bahahascf_9", 'te_fraction'])
datHighImpactVsTENorm <- cbind(datHighImpactVsTE[, c('chr', 'start', 'end', 'te_fraction','dna_te_fraction', 'ltr_fraction')], datHighImpactVsTE[, 8:ncol(datHighImpactVsTE)] / datHighImpactVsTE$cds_fraction);

#hist(datHighImpactVsTE$cds_fraction);
hist(log10(datHighImpactVsTE$te_fraction) );
hist(log10(datHighImpactVsTENorm$af0_20+1e-5) )
fnWilcox <- function(sTEType, sAF) {
  arrNoHighImpact <- datHighImpactVsTENorm[datHighImpactVsTENorm[,sAF]==0, sTEType]
  arrWithHighImpact <- datHighImpactVsTENorm[datHighImpactVsTENorm[,sAF]>0, sTEType ]
  arrNoHighImpact <- arrNoHighImpact[!is.na(arrNoHighImpact)]
  arrWithHighImpact <- arrWithHighImpact[!is.na(arrWithHighImpact)]
  cat("median ", sTEType ," in regions with no high impact variants: ", median(arrNoHighImpact), "\n");
  cat("median ", sTEType ," in regions with high impact variants: ", median(arrWithHighImpact), "\n");
  print(wilcox.test(arrNoHighImpact, arrWithHighImpact))  
  cat("mean ", sTEType ," in regions with no high impact variants: ", mean(arrNoHighImpact), "\n");
  cat("mean ", sTEType ," in regions with high impact variants: ", mean(arrWithHighImpact), "\n");
  print(t.test(arrNoHighImpact, arrWithHighImpact))  
  
}

datHighImpactVsTEFullModel <- NULL;
for(sAF in colnames(datHighImpactVsTE)[8:ncol(datHighImpactVsTE)] ) {
  cat("===========allele freq range ", sAF, " Excluded wins with zero high impact var: ", bExcludeZeros, " ==================\n");
  
  #plot(datHighImpactVsTENorm[,sAF] ~ log10(datHighImpactVsTENorm[, 'te_fraction']), pch=16, cex=0.3);  
  if (!bExcludeZeros) {
    fnWilcox('te_fraction', sAF);
    fnWilcox('dna_te_fraction', sAF);
    fnWilcox('ltr_fraction', sAF);
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm;
    datHighImpactVsTENorm[,sAF] <- datHighImpactVsTENorm[,sAF] + 1e-5; #add small value for log transformation.
  } else {
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm[datHighImpactVsTENorm[,sAF]>0 ,];
  }
  
  
  print(oSum1<-summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'te_fraction']) )));
  print(oSum2<-summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'dna_te_fraction']) )));
  print(oSum3<-summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'ltr_fraction']) )));

  plot(log10(datHighImpactVsTENorm2[, 'te_fraction']) , log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#3a9ab299", main=sAF, xlab="log10(TE fraction)", ylab="log10(Normalized high impact var. count)");
  points(log10(datHighImpactVsTENorm2[, 'dna_te_fraction']), log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#e1bc2199");
  points(log10(datHighImpactVsTENorm2[, 'ltr_fraction']), log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1,col="#f11b0099");
  
  abline(coef = oSum1$coefficients, col="#3a9ab2ff", lwd=3)
  abline(coef = oSum2$coefficients, col="#e1bc21ff", lwd=3)
  abline(coef = oSum3$coefficients, col="#f11b00ff", lwd=3)
  
  datHighImpactVsTEFullModel <- rbind(datHighImpactVsTEFullModel, data.frame(rho=datHighImpactVsTENorm2[,sAF], 
                                                                             AF = sAF, 
                                                                             te_fraction=datHighImpactVsTENorm2[, 'te_fraction'], 
                                                                             dna_te_fraction=datHighImpactVsTENorm2[, 'dna_te_fraction'],
                                                                             ltr_fraction = datHighImpactVsTENorm2[, 'ltr_fraction']
                                                                             ) );

}

summary(lm(log10(datHighImpactVsTEFullModel$rho) ~ log10(datHighImpactVsTEFullModel[, 'te_fraction']) + datHighImpactVsTEFullModel[, 'AF'] ))
summary(lm(log10(datHighImpactVsTEFullModel$rho) ~ log10(datHighImpactVsTEFullModel[, 'dna_te_fraction']) + datHighImpactVsTEFullModel[, 'AF'] ))
summary(lm(log10(datHighImpactVsTEFullModel$rho) ~ log10(datHighImpactVsTEFullModel[, 'ltr_fraction']) + datHighImpactVsTEFullModel[, 'AF'] ))

datChrLen <- aggregate(datHighImpactVsTENorm[, c('end')], by = list(chr=datHighImpactVsTENorm$chr), FUN=max)
colnames(datChrLen) <- c('chr', 'chrlen');
datHighImpactVsTENorm$winmid <- (datHighImpactVsTENorm$start + datHighImpactVsTENorm$end)/2
datHighImpactVsTENorm <- merge(datHighImpactVsTENorm, datChrLen, by="chr", all.x =T, all.y=F);
datHighImpactVsTENorm$winrelativepos <- datHighImpactVsTENorm$winmid / datHighImpactVsTENorm$chrlen;
datHighImpactVsTENorm$winrelativepos[datHighImpactVsTENorm$winrelativepos>0.5] <- 1 - datHighImpactVsTENorm$winrelativepos[datHighImpactVsTENorm$winrelativepos>0.5];

datHighImpactVsTEFullModel <- NULL;
for(sAF in colnames(datHighImpactVsTE)[8:ncol(datHighImpactVsTE)] ) {
  cat("===========Pos. allele freq range ", sAF, " Excluded wins with zero high impact var: ", bExcludeZeros, " ==================\n");
  
  #plot(datHighImpactVsTENorm[,sAF] ~ log10(datHighImpactVsTENorm[, 'te_fraction']), pch=16, cex=0.3);  
  if (!bExcludeZeros) {
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm;
    datHighImpactVsTENorm[,sAF] <- datHighImpactVsTENorm[,sAF] + 1e-5; #add small value for log transformation.
  } else {
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm[datHighImpactVsTENorm[,sAF]>0 ,];
  }
  
  
  print(oSum1<-summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ datHighImpactVsTENorm2[, 'winrelativepos']) ));

  plot( datHighImpactVsTENorm2[, 'winrelativepos'] , log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#3a9ab299", main=sAF, xlab="Relative distance to nearest telomere", ylab="log10(Normalized high impact var. count)");
  
  abline(coef = oSum1$coefficients, col="#3a9ab2ff", lwd=3)
  datHighImpactVsTEFullModel <- rbind(datHighImpactVsTEFullModel, data.frame(rho=datHighImpactVsTENorm2[,sAF], 
                                                                             AF = sAF, 
                                                                             winrelativepos=datHighImpactVsTENorm2[, 'winrelativepos']
  ) );
  
}

print(oSum1<-summary(lm(log10(datHighImpactVsTENorm[,'te_fraction']) ~ datHighImpactVsTENorm[, 'winrelativepos']) ));

plot( datHighImpactVsTENorm[, 'winrelativepos'] , log10(datHighImpactVsTENorm[,'te_fraction']), pch=16, cex=1, col="#3a9ab299", main="TE fraction", xlab="Relative distance to nearest telomere", ylab="log10(TE fraction)");

abline(coef = oSum1$coefficients, col="#3a9ab2ff", lwd=3)

print(oSum1<-summary(lm(log10(datHighImpactVsTEFullModel[,'rho']) ~ datHighImpactVsTEFullModel[, 'winrelativepos'] + datHighImpactVsTEFullModel[, 'AF'] ) ));


dev.off();