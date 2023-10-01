setwd("/public3/group_crf/home/cuirf/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/rho");
source("asym.MK.ray.R");

arrImpact <- c("missense_variant");
lsPlotData <- list();

sOut <- "asymtotoic_mk_nonsyn.txt"
arrBreaks <- c(seq(0.0, 0.80, 0.05) , 0.99, 1);
arrAlleleFreq <- arrBreaks[2:(length(arrBreaks)-1)]
AF_low <- 0.08
AF_high <- 0.92
#dat <- read.table("consurf_af_sv.rec.bed", header=F, sep="\t")
dat <- read.table("NONSYN.polarized.out.txt", header=F, sep="\t")

dat <- dat[ , c(1,2,2, 3:ncol(dat) ) ]
dat[,2] <- dat[,2] - 1;
dat <- dat[complete.cases(dat), ];
dat$AF <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; return(p) })
dat <- dat[dat$AF>0, ];
dat$AFBin <- cut(dat$AF, breaks = arrBreaks);
dat <- dat[complete.cases(dat),];

sink(file = sOut)

for(sImpact in arrImpact) {
  #plot asymtotic MK-like :
  dN <- sum(as.integer(grepl( sImpact, dat$V6 ) & dat$AF == 1));
  dS <- sum(as.integer(grepl(  "synonymous_variant", dat$V6 ) & dat$AF == 1));
  pN <- sum(as.integer(grepl( sImpact, dat$V6 ) & dat$AF < 1));
  pS <- sum(as.integer(grepl(  "synonymous_variant", dat$V6 ) & dat$AF < 1));
  
  arrPN <- c();
  arrPS <- c();
  arrAFBins <- as.character(levels(dat$AFBin));
  
  for (sAFBin in arrAFBins[1:(length(arrAFBins)-1)]) {
    arrPN <- c(arrPN, sum(as.integer(grepl( sImpact, dat$V6 ) & dat$AFBin == sAFBin )));
    arrPS <- c(arrPS, sum(as.integer(grepl(  "synonymous_variant", dat$V6 ) & dat$AFBin == sAFBin )));
  }
  

  alpha <- 1 - ( (dS * pN) / (dN * pS) );
  oConTable <- matrix(c(dN , dS, pN , pS ), nrow = 2);
  oTest <- fisher.test(oConTable);
  cat(sImpact, ": ===============\n")
  cat("dN\tdS\tpN\tpS\n");
  cat(dN, "\t", dS, "\t", pN, "\t", pS, "\n");
  cat("alpha global = ", alpha);
  print(oTest)
  
  arrAsymMK <- asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=AF_low, xhigh=AF_high, output='table')
  print(arrAsymMK)
  
  n10percIdx <- which(arrAlleleFreq == 0.1);
  alpha_10perc <- 1 - ( (dS * arrPN[n10percIdx] ) / (dN * arrPS[n10percIdx] ) );
  alpha_asymptotic <- arrAsymMK['alpha_asymptotic'];
  
  alphaDiff <- alpha_asymptotic - alpha_10perc;
  
  cat("alpha_10%\talpha_original\talpha_asymptotic\talphaDiff\n");
  cat(alpha_10perc, "\t",  as.numeric(arrAsymMK['alpha_original']) , "\t", as.numeric(alpha_asymptotic), "\t", as.numeric(alphaDiff), "\n");

  lsPlotData[[sImpact]] <- list();
  lsPlotData[[sImpact]]$d0 <- dS;
  lsPlotData[[sImpact]]$d <- dN;
  lsPlotData[[sImpact]]$xlow <- AF_low;
  lsPlotData[[sImpact]]$xhigh <- AF_high;
  lsPlotData[[sImpact]]$df <- data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL);
  lsPlotData[[sImpact]]$true_alpha <- NA;
  lsPlotData[[sImpact]]$arrYLims <- c(-2, 1);
  
}

sink();


#lsPlotData[[1]]$col <- '#f21a00';
#lsPlotData[[2]]$col <- '#ebcc2a';
lsPlotData[[1]]$col <- '#3b9ab2';


pdf(file=paste(sOut, ".pdf",sep="") , width=6, height=4)
#asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.08, xhigh=0.9, arrYLims=c(-1, 0.5), output='default')

asymptoticMKOverlay(lsPlotData);
#abline(h=alpha_10perc , col='blue');
#abline(h= as.numeric(arrAsymMK['alpha_original'])  , col='green');
#abline(h= 0  , col='black');




##plot mean consurf scores
oF <- gzfile("/public3/group_crf/home/cuirf/bahaha_assembly/consurf/bahaha/refgenome/consurf_scores_flat.txt.gz")
datConsurf <- read.table(oF, sep="\t", header=F)
datConsurf$coord <- paste(datConsurf$V1, datConsurf$V2, sep=":");

datMissense <- dat[grepl( "missense_variant", dat$V6 ), ];
datMissense$coord <- paste(datMissense$V1, datMissense$V2.1, sep=":");
colnames(datConsurf)[3] <- "ConsurfScore"
datMissense <- merge(datMissense,datConsurf[,c(3,5)], by="coord" )

arrBreaks <- c(seq(0.0, 0.8, 0.2) , 0.99, 1);


datMissense$AFBin <- cut(datMissense$AF, breaks = arrBreaks);


library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sumMoreThan8 = sum(x[[col]][x[[col]]>=8], na.rm=TRUE),
      sumEqual1 = sum(x[[col]][x[[col]]==1], na.rm=TRUE),
      sumEqual2 = sum(x[[col]][x[[col]]==2], na.rm=TRUE),
      sumEqual3 = sum(x[[col]][x[[col]]==3], na.rm=TRUE),
      sumEqual4 = sum(x[[col]][x[[col]]==4], na.rm=TRUE),
      sumEqual5 = sum(x[[col]][x[[col]]==5], na.rm=TRUE),
      sumEqual6 = sum(x[[col]][x[[col]]==6], na.rm=TRUE),
      sumEqual7 = sum(x[[col]][x[[col]]==7], na.rm=TRUE),
      sumEqual8 = sum(x[[col]][x[[col]]==8], na.rm=TRUE),
      sumEqual9 = sum(x[[col]][x[[col]]==9], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]][!is.na(x[[col]])] )))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


dsum <- data_summary(datMissense, 'ConsurfScore'  , c('AFBin') );
dsum
summary(lm(datMissense$ConsurfScore ~ datMissense$AF ));

arrYLims <- c(4,6)
p <- ggplot(dsum, aes(x=AFBin, y=ConsurfScore, fill="#08bdbd")) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=ConsurfScore-sem, ymax=ConsurfScore+sem), width=.2,
                position=position_dodge(.9)) + coord_cartesian(ylim = arrYLims )

p + scale_fill_manual(values=c( "#08bdbd")) + theme_classic()

dev.off();
