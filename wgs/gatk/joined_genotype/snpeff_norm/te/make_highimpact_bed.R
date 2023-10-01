setwd("~/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/te")
library(vioplot)
library(stringr)
library(ggplot2)
library(wesanderson)
arrBreaks <- c(seq(0.0, 0.8, 0.2) , 0.99, 1);

dat <- read.table("../rho/consurf_af_sv.rec.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
dat <- dat[complete.cases(dat),];

datHigh <- dat[dat$V10 == 'HIGH',];

datHigh$AF <- apply((datHigh[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; return(p) })

datHigh <- datHigh[datHigh$AF > 0, ]
datHigh$AFBin <- cut(datHigh$AF, breaks = arrBreaks);
for(sAFBin in levels(datHigh$AFBin)) {
  sRange <- str_replace_all(str_replace_all(sAFBin, '[\\(\\])]', ''), ',','-');
  sOutFile <- paste0("highimpact.af",sRange,".bed");
  datWrite <- datHigh[datHigh$AFBin == sAFBin, c(1,2,3)];
  write.table(datWrite, file=sOutFile, col.names = F, row.names = F, quote = F, sep="\t");
}