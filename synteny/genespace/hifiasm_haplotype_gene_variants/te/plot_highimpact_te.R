setwd("~/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/te")

pdf(file="te_correlate.pdf", width = 5, height=5)
datHighImpactVsTE <- read.table("cov.perwin.bed", header=F, sep="\t")
colnames(datHighImpactVsTE) <- c('chr', 'start', 'end', 'te_fraction' , 'dna_te_fraction', 'ltr_fraction', 'cds_fraction', 'allhaps', 'F0', 'F1', 'M0', 'M1')

# nMinCDSFraction <- quantile(datHighImpactVsTE$cds_fraction, c(0.025))
# datHighImpactVsTE <- datHighImpactVsTE[datHighImpactVsTE$cds_fraction >=nMinCDSFraction,];

bExcludeZeros <- T;

plot(datHighImpactVsTE[datHighImpactVsTE$chr=="bahahascf_9", 'start'],datHighImpactVsTE[datHighImpactVsTE$chr=="bahahascf_9", 'te_fraction'])
datHighImpactVsTENorm <- cbind(datHighImpactVsTE[, c('chr', 'start', 'end', 'te_fraction','dna_te_fraction', 'ltr_fraction')], datHighImpactVsTE[, 8:ncol(datHighImpactVsTE)] / datHighImpactVsTE$cds_fraction);

#hist(datHighImpactVsTE$cds_fraction);
hist(log10(datHighImpactVsTE$te_fraction) );
hist(log10(datHighImpactVsTENorm$F0) )
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

for(sAF in colnames(datHighImpactVsTE)[8:ncol(datHighImpactVsTE)] ) {
  cat("===========Haplotype: ", sAF, " Excluded wins with zero high impact var: ", bExcludeZeros, " ==================\n");
  
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
  
  
  print(oSum1 <- summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'te_fraction']) )));
  print(oSum2 <- summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'dna_te_fraction']) )));
  print(oSum3 <- summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ log10(datHighImpactVsTENorm2[, 'ltr_fraction']) )));
  
  plot(log10(datHighImpactVsTENorm2[, 'te_fraction']) , log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#3a9ab299", main=sAF, xlab="log10(TE fraction)", ylab="log10(Normalized high impact var. count)");
  points(log10(datHighImpactVsTENorm2[, 'dna_te_fraction']), log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#e1bc2199");
  points(log10(datHighImpactVsTENorm2[, 'ltr_fraction']), log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1,col="#f11b0099");

  abline(coef = oSum1$coefficients, col="#3a9ab2ff", lwd=3)
  abline(coef = oSum2$coefficients, col="#e1bc21ff", lwd=3)
  abline(coef = oSum3$coefficients, col="#f11b00ff", lwd=3)
  
}

datChrLen <- aggregate(datHighImpactVsTENorm[, c('end')], by = list(chr=datHighImpactVsTENorm$chr), FUN=max)
colnames(datChrLen) <- c('chr', 'chrlen');
datHighImpactVsTENorm$winmid <- (datHighImpactVsTENorm$start + datHighImpactVsTENorm$end)/2
datHighImpactVsTENorm <- merge(datHighImpactVsTENorm, datChrLen, by="chr", all.x =T, all.y=F);
datHighImpactVsTENorm$winrelativepos <- datHighImpactVsTENorm$winmid / datHighImpactVsTENorm$chrlen;
datHighImpactVsTENorm$winrelativepos[datHighImpactVsTENorm$winrelativepos>0.5] <- 1 - datHighImpactVsTENorm$winrelativepos[datHighImpactVsTENorm$winrelativepos>0.5];


for(sAF in colnames(datHighImpactVsTE)[8:ncol(datHighImpactVsTE)] ) {
  cat("===========Pos. haplotype ", sAF, " Excluded wins with zero high impact var: ", bExcludeZeros, " ==================\n");
  
  #plot(datHighImpactVsTENorm[,sAF] ~ log10(datHighImpactVsTENorm[, 'te_fraction']), pch=16, cex=0.3);  
  if (!bExcludeZeros) {
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm;
    datHighImpactVsTENorm[,sAF] <- datHighImpactVsTENorm[,sAF] + 1e-5; #add small value for log transformation.
  } else {
    datHighImpactVsTENorm2 <- datHighImpactVsTENorm[datHighImpactVsTENorm[,sAF]>0 ,];
  }
  
  
  print(oSum1 <-summary(lm(log10(datHighImpactVsTENorm2[,sAF]) ~ datHighImpactVsTENorm2[, 'winrelativepos']) ));
  plot( datHighImpactVsTENorm2[, 'winrelativepos'] , log10(datHighImpactVsTENorm2[,sAF]), pch=16, cex=1, col="#3a9ab299", main=sAF, xlab="Relative distance to nearest telomere", ylab="log10(Normalized high impact var. count)");
  
  abline(coef = oSum1$coefficients, col="#3a9ab2ff", lwd=3)
}

dev.off();
