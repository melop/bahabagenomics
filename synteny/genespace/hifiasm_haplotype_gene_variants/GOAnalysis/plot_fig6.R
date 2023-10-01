setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/GOanalysis");
library(chromPlot)
library(RColorBrewer)
library(wesanderson)

arrChrs <- 1:24;
nMinWinPlotLen <- 3e5;
sRefSp <- "BahahaTaipingensis"

datRNACoord <- read.table("mRNA.coord.txt", header=T, sep="\t", quote="")
datOrthID2RNAID <- read.table("../../hifiasm_bahaha_haplotypes/synorthos.txt", header=T, sep="\t")
datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)

datP <- datPresence[,5:ncol(datPresence)]
datH <- datHighImpact[,5:ncol(datHighImpact)]

nSamples <- ncol(datP)
arrSamples <- colnames(datP)
datRet <- datH

for(nIdx in 1:ncol(datP) ) {
  datRet[, nIdx] <- as.integer(datP[ , nIdx] >=0.5 & datH[, nIdx]==0)
}


datOnlyCompHighImpact <- cbind(datPresence[,1:3], datRet);
datOnlyCompHighImpact$UniqID <- paste(datOnlyCompHighImpact$pgChr, datOnlyCompHighImpact$pgOrd, datOnlyCompHighImpact$pgID, sep="_")
datOrthID2RNAID$UniqID <- paste(datOrthID2RNAID$pgChr, datOrthID2RNAID$pgOrd, datOrthID2RNAID$pgID, sep="_")
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datOrthID2RNAID[, c("UniqID",  sRefSp)], by="UniqID", all.x = T, all.y=F);

arrNonRepresentativeInDupArray <- datOnlyCompHighImpact$BahahaTaipingensis[! (datOnlyCompHighImpact$BahahaTaipingensis %in% datRNACoord$mRNAID) ]

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datRNACoord, by.x = sRefSp, by.y="mRNAID");

summary(lm(datOnlyCompHighImpact$BahahaTaipingensisF0 ~ datOnlyCompHighImpact$BahahaTaipingensisF1))

arrHapNames <- colnames(datOnlyCompHighImpact)[6:9]
datHyperGeoTests <- NULL;
datHyperGeoTestsPerChr  <- NULL;
datPerHapPerChrDel <- NULL;

for(nHap1 in 1:(length(arrHapNames)) ) {
  sHap1 <- arrHapNames[nHap1];
  for(nChr in arrChrs) {
    
    nCount <- sum(as.integer(datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr, sHap1]  == 0 ))
    datPerHapPerChrDel <- rbind(datPerHapPerChrDel, data.frame(hap1=sHap1, chr=nChr, DelCount=nCount));
    
  }
}

for(nHap1 in 1:(length(arrHapNames)-1) ) {
  sHap1 <- arrHapNames[nHap1];
  
  
  
  for(nHap2 in (nHap1+1):length(arrHapNames) ) {
    sHap2 <- arrHapNames[nHap2];

    cat("\nbetween ", sHap1, " ", sHap2 ,"\n\n");
    
    nBothOK <- sum(as.integer(datOnlyCompHighImpact[ , sHap1]  == 1 & datOnlyCompHighImpact[ , sHap2] == 1))
    nHap2Del <- sum(as.integer(datOnlyCompHighImpact[ , sHap1]  == 1 & datOnlyCompHighImpact[ , sHap2] == 0))
    nHap1Del <- sum(as.integer(datOnlyCompHighImpact[ , sHap1]  == 0 & datOnlyCompHighImpact[ , sHap2] == 1))
    nBothHapDel <- sum(as.integer(datOnlyCompHighImpact[ , sHap1]  == 0 & datOnlyCompHighImpact[ , sHap2] == 0))
    datHyperGeoTests <- rbind(datHyperGeoTests, data.frame(hap1=sHap1, hap2=sHap2, BothOK=nBothOK, hap1onlydel=nHap1Del, hap2onlydel=nHap2Del, bothhapdel=nBothHapDel,
    hyperP=phyper(nBothHapDel, nHap2Del + nBothHapDel, nHap1Del + nBothOK, nHap1Del + nBothHapDel, lower.tail = FALSE)));
    
    for(nChr in arrChrs) {
      cat("\nbetween ", sHap1, " ", sHap2 ,"\n\n");
      
      nBothOK <- sum(as.integer(datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr, sHap1]  == 1 & datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap2] == 1))
      nHap2Del <- sum(as.integer(datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap1]  == 1 & datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap2] == 0))
      nHap1Del <- sum(as.integer(datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap1]  == 0 & datOnlyCompHighImpact[ datOnlyCompHighImpact$pgChr ==  nChr , sHap2] == 1))
      nBothHapDel <- sum(as.integer(datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap1]  == 0 & datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr ==  nChr , sHap2] == 0))
      datHyperGeoTestsPerChr <- rbind(datHyperGeoTestsPerChr, data.frame(hap1=sHap1, hap2=sHap2, chr=nChr, BothOK=nBothOK, hap1onlydel=nHap1Del, hap2onlydel=nHap2Del, bothhapdel=nBothHapDel,
                                                             hyperP=phyper(nBothHapDel, nHap2Del + nBothHapDel, nHap1Del + nBothOK, nHap1Del + nBothHapDel, lower.tail = FALSE)));
      
    }
  }
}

datHyperGeoTests$hap1overlappercent <- datHyperGeoTests$bothhapdel / (datHyperGeoTests$hap1onlydel + datHyperGeoTests$bothhapdel)
datHyperGeoTests$hap2overlappercent <- datHyperGeoTests$bothhapdel / (datHyperGeoTests$hap2onlydel + datHyperGeoTests$bothhapdel)
write.table(datHyperGeoTests, file="hypergeometric_test_overlap_hap_delvars.txt", col.names = T, row.names = F, quote = F, sep="\t");

datHyperGeoTestsPerChr$hap1overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap1onlydel + datHyperGeoTestsPerChr$bothhapdel)
datHyperGeoTestsPerChr$hap2overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap2onlydel + datHyperGeoTestsPerChr$bothhapdel)
write.table(datHyperGeoTestsPerChr, file="hypergeometric_test_overlap_hap_delvars_perchr.txt", col.names = T, row.names = F, quote = F, sep="\t");

#draw variation between haplotypes
# dat2Bind <- cbind( datHyperGeoTestsPerChr[ , c(3,4)] ,'BothOK' );
# names(dat2Bind)[2:3] <- c("count", "type");
# datHyperGeoTestsPerChrPlot <- dat2Bind; 
dat2Bind <- cbind( datPerHapPerChrDel[ , c(2,3)] ,'del_per_hap' );
names(dat2Bind)[2:3] <- c("count", "type");
datHyperGeoTestsPerChrPlot <- dat2Bind; 
dat2Bind <- cbind( datHyperGeoTestsPerChr[ , c(3,5)] ,'private_to_hap' );
names(dat2Bind)[2:3] <- c("count", "type");
datHyperGeoTestsPerChrPlot <- rbind(datHyperGeoTestsPerChrPlot, dat2Bind); 
dat2Bind <-cbind( datHyperGeoTestsPerChr[ , c(3,6)] ,'private_to_hap' ) ;
names(dat2Bind)[2:3] <- c("count", "type");
datHyperGeoTestsPerChrPlot <- rbind(datHyperGeoTestsPerChrPlot, dat2Bind); 
dat2Bind <- cbind( datHyperGeoTestsPerChr[ , c(3,7)] ,'shared_between_haps' );
names(dat2Bind)[2:3] <- c("count", "type");
datHyperGeoTestsPerChrPlot <- rbind(datHyperGeoTestsPerChrPlot, dat2Bind); 

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  colnames(datac)[colnames(datac)=='mean'] <- measurevar;
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


datHapDelPerChrPlot <- summarySE(data=datHyperGeoTestsPerChrPlot, measurevar="count", groupvars=c("chr","type"))
ggplot(datHapDelPerChrPlot, aes(fill=type, y=count, x=chr)) + 
  geom_bar(position='dodge',  stat="identity") +
  geom_errorbar(aes(ymin=count-sd,
                    ymax=count+sd),
               size=0.2, width=0.7, linetype=1, position=position_dodge(.9)) +
  theme_classic() 
ggsave("del_sd_per_chr.pdf", width=12, height=2.5)

sChrFile <- "chromosomes.tsv";


#chrom plot:
datChr <- read.table(sChrFile, sep="\t", header = F)
datChr$Name <- datChr$V1
colnames(datChr) <- c("Chrom", "Start", "End", "Name");
datChr$Chrom <- gsub("bahahascf_","", datChr$Chrom)
datChr$Start <- datChr$End - 10000;
datChr$Name <- 'telomere';

datChr2 <- datChr;
datChr2$Start <- 0;
datChr2$End <- 10000;

# datChr3 <- datChr;
# datChr3$Start <- as.integer(datChr3$End/2) -1000;
# datChr3$End <- as.integer(datChr3$End/2);
# datChr3$Name <- 'centromere'

datChr <- rbind(datChr, datChr2);
datChr <- datChr[order(datChr$Chrom, datChr$Start), ]

datFakeMarkers <- datChr[, 1:3]
datFakeMarkers$Group <- 'telomere'
datFakeMarkers$Colors <- '#FFFFFF00'

#

datPerHapPerChrDel <- merge(datPerHapPerChrDel, datChr, by.x="chr", by.y="Chrom", all.x = T, all.y =F);
datPerHapPerChrDelForTest <- datPerHapPerChrDel[datPerHapPerChrDel$Start>0, ]
summary(lm(datPerHapPerChrDelForTest$DelCount ~ datPerHapPerChrDelForTest$End))
arrGeneCountsPerChr <- table(datPresence$pgChr)[1:24]
datGeneCountsPerChr <- data.frame(chr=names(arrGeneCountsPerChr), genecount=arrGeneCountsPerChr )
datPerHapPerChrDelForTest <- merge(datPerHapPerChrDelForTest, datGeneCountsPerChr, by.x="chr", by.y="chr", all.x = T, all.y =F);
summary(lm(datPerHapPerChrDelForTest$DelCount ~ datPerHapPerChrDelForTest$genecount.Freq ))

datDelCountFoldDiffHaps <- aggregate(datPerHapPerChrDel$DelCount ,by = list(chr=datPerHapPerChrDel$chr) , FUN=function(x) {max(x)/min(x)})

arrHaplotypeIDs <- colnames(datOnlyCompHighImpact)[1+grep(sRefSp, colnames(datOnlyCompHighImpact)[2:ncol(datOnlyCompHighImpact)])]
datROHSegments <- NULL;
arrColors <- (as.character(wes_palette("Zissou1", length(arrHaplotypeIDs), type = "continuous")));
arrColors <- paste0(arrColors, 'aa');
nColor <- 1;
for (sF in arrHaplotypeIDs ) {
  
  datThis <- datOnlyCompHighImpact[datOnlyCompHighImpact[,sF] == 0, c("chr", "start", "end") ]
  datThis$Group <- sF;
  datThis$Colors <- arrColors[nColor]
  nColor <- nColor + 1;
  datROHSegments <- rbind(datROHSegments , datThis);
}


#CHECK IF DISTRIBUTION OF VARIANTS ARE CLUSTERED
datROHSegments$mid <-  ( datROHSegments$end + datROHSegments$start )/2;
arrDist <- c();
datVarCounts <- NULL;
for(sInd in unique(datROHSegments$Group)) {
  for(sChr in unique(datROHSegments$chr)) {
    arrMids <- datROHSegments[datROHSegments$chr==sChr & datROHSegments$Group == sInd, 'mid'];
    datVarCounts <- rbind(datVarCounts, data.frame(ind=sInd, chr=sChr, count=length(arrMids) ) );
    for(i in 1:length(arrMids)) {
      for(j in 1:length(arrMids)) {
        if (i>=j) {next;}
        arrDist <- c(arrDist, abs(arrMids[i] - arrMids[j]) );
      }
    }
  }
}


#hist(arrDist);


lsDistSim <- list();
datRNACoord$mid <- (datRNACoord$start + datRNACoord$end)/2;
for(nRep in 1:1000) {
  arrDistSim <- c();
  cat("Rep ", nRep, "\n");
  for(sInd in unique(datVarCounts$ind)) {
    for(sChr in unique(datVarCounts$chr)) {
      nChrLen <- datChr[datChr$Name == sChr, 'End']
      if (length(nChrLen) == 0) {next;}
      nDelCountPerChr <- datVarCounts[datVarCounts$ind==sInd & datVarCounts$chr==sChr, 'count' ];
      if (length(nDelCountPerChr) == 0) {next;}
      
      arrMidChr <- datRNACoord[datRNACoord$chr==sChr, 'mid']
      #arrMids <- as.integer(runif(nDelCountPerChr, 1, nChrLen ));
      arrMids <- sample(x = arrMidChr , size = nDelCountPerChr, replace = F);
      for(i in 1:length(arrMids)) {
        for(j in 1:length(arrMids)) {
          if (i>=j) {next;}
          arrDistSim <- c(arrDistSim, abs(arrMids[i] - arrMids[j]) );
        }
      }
    }
  }
  #arrDistSimAvg <- c(arrDistSimAvg, mean(arrDistSim[arrDistSim<5e5]) );
  lsDistSim[[nRep]] <-  arrDistSim;
}

datDistP <- NULL;
for(nDistCutoffs in seq(1e5, 40e6, 1e5) ) {
  nObs <- sd(arrDist[arrDist<nDistCutoffs]);
  arrSim <- unlist(lapply(lsDistSim, function(x) {sd(x[x<nDistCutoffs])} ))
  arrSimQuantiles <- quantile(arrSim, c(0.005, 0.025, 0.5, 0.975, 0.995) )
  datQuantiles <- data.frame(q005 = arrSimQuantiles[1], q025= arrSimQuantiles[2], q5= arrSimQuantiles[3], q975= arrSimQuantiles[4], q995= arrSimQuantiles[5])
  nLowCount <- length(arrSim[arrSim<nObs]);
  if (nLowCount==0) {
    nLowCount <- 0.5
  }
  if (nLowCount==length(arrSim)) {
    nLowCount <- length(arrSim) - 0.5
  }
  
  datDistPThis<- data.frame(distcutoff=nDistCutoffs, obs=nObs, p=nLowCount/length(arrSim) );
  datDistPThis <- cbind(datDistPThis, datQuantiles);
  datDistP <- rbind(datDistP, datDistPThis);
}

datDistP$p_twotail <- datDistP$p
datDistP$p_twotail[datDistP$p<0.5] <- datDistP$p[datDistP$p<0.5] * 2;
datDistP$p_twotail[datDistP$p>=0.5] <- (1- datDistP$p[datDistP$p>=0.5]) * 2;

plot(datDistP$distcutoff, -log10(datDistP$p_twotail), type='l' )

arrSim <- unlist(lapply(lsDistSim, function(x) {sd(x)} ))
nObs <- sd(arrDist);

mean(arrDistNoLong);
quantile(arrDistSimAvg, c(0.005, 0.025, 0.975, 0.995));
pdf(file="del_dist_sd_sim_vs_obs.pdf", width = 4,height = 6)
hist(arrSim, xlim = c(6400000,8000000 ), col = "#9EBE91aa" );
abline(v=nObs, col="red")
length(arrSim[arrSim>nObs]) / length(arrSim);
dev.off();
#simulate:


#

datROHSegments$len <- datROHSegments$end - datROHSegments$start + 1
arrTooSmallWin <- datROHSegments$len<nMinWinPlotLen;
datROHSegments$start[arrTooSmallWin] <- datROHSegments$start[arrTooSmallWin] - (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2
datROHSegments$end[arrTooSmallWin] <- datROHSegments$end[arrTooSmallWin] + (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2

datROHSegmentsForPlot <- data.frame(Chrom=gsub("bahahascf_","", datROHSegments$chr), Start=datROHSegments$start, End=datROHSegments$end, Group=factor(x = datROHSegments$Group, levels = unique(datROHSegments$Group)), Colors=factor(datROHSegments$Colors, levels = arrColors) )

datROHSegmentsForPlot <- datROHSegmentsForPlot[as.integer(datROHSegmentsForPlot$Chrom)<=24, ]


datOnlyCompHighImpact$sumDel <- rowSums(datOnlyCompHighImpact[, 1+grep(sRefSp, colnames(datOnlyCompHighImpact)[2:ncol(datOnlyCompHighImpact)]) ])
datGeneLabels <- datOnlyCompHighImpact[datOnlyCompHighImpact$sumDel<3, c('chr', 'start', 'end', 'GeneSymbol', 'sumDel')]
datGeneLabels <- datGeneLabels[(!grepl('unknown', datGeneLabels$GeneSymbol, fixed = T)) & (!grepl('loc[0-9]*',datGeneLabels$GeneSymbol )) & (!grepl('si:',datGeneLabels$GeneSymbol )), ];
datGeneLabelsPlot <- data.frame(Chrom=gsub("bahahascf_","", datGeneLabels$chr), Start=datGeneLabels$start, End=datGeneLabels$end, ID=datGeneLabels$GeneSymbol );
pdf(file = "fig6.highimpact.sv.del.pdf",width=10/1.5,height=8/1.5);

datROHSegmentsForPlot <- rbind(datROHSegmentsForPlot, datFakeMarkers)
chromPlot(gaps=datChr, bands=datROHSegmentsForPlot, figCols=12, stat=datGeneLabelsPlot , 
          statCol="Value",  
          noHist = T,
          statName ="genes", statTyp="n", statSumm="none" )




dev.off();


arrPerc <- (datHyperGeoTests$hap1onlydel + datHyperGeoTests$bothhapdel) / (datHyperGeoTests$BothOK + datHyperGeoTests$hap1onlydel + datHyperGeoTests$hap2onlydel + datHyperGeoTests$bothhapdel)
arrPerc <- c(arrPerc, (datHyperGeoTests$hap2onlydel + datHyperGeoTests$bothhapdel) / (datHyperGeoTests$BothOK + datHyperGeoTests$hap1onlydel + datHyperGeoTests$hap2onlydel + datHyperGeoTests$bothhapdel))
mean(arrPerc)
quantile(arrPerc, c(0.025, 0.975) );
