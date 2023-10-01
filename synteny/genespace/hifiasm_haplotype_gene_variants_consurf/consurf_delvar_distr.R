setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf");
#analyze the distribution of high consurf scores on 4 haplotypes
library(chromPlot)
library(RColorBrewer)
library(ggplot2)
library(wesanderson)


datConsurfRec <- read.table("consurf_af.rec.bed", header=F, sep="\t")
datConsurfRec$chr <- as.integer(gsub("bahahascf_", "", datConsurfRec$V1))

datConsurfRec <- datConsurfRec[ rowSums(datConsurfRec[,6:9])<4, ];
datConsurfRec <- datConsurfRec[complete.cases(datConsurfRec[,6:9]),];#only keep segregating variants

datConsurfRec1 <- datConsurfRec[datConsurfRec$V5==1, ];
datConsurfRec9 <- datConsurfRec[datConsurfRec$V5>=9, ];

mean(datConsurfRec1$V13)
mean(datConsurfRec9$V13)

wilcox.test(datConsurfRec1$V13,datConsurfRec9$V13 )

# 
# sum(datConsurfRec9$V6 == 1, na.rm=T)
# sum(datConsurfRec9$V7 == 1, na.rm=T)
# sum( (datConsurfRec9$V6 + datConsurfRec9$V7) == 2, na.rm = T)

datRNACoord <- read.table("../hifiasm_haplotype_gene_variants/GOanalysis/mRNA.coord.txt", header=T, sep="\t", quote="")
datConsurfRec9 <- merge(datConsurfRec9, datRNACoord, by.x = 'V4', by.y="mRNAID");


arrChrs <- 1:24
colnames(datConsurfRec9)[6:9] <- c('F0', 'F1', 'M0', 'M1');
arrHapNames <- colnames(datConsurfRec9)[6:9]
datHyperGeoTests <- NULL;
datHyperGeoTestsPerChr  <- NULL;
datPerHapPerChrDel <- NULL;

for(nHap1 in 1:(length(arrHapNames)) ) {
  sHap1 <- arrHapNames[nHap1];
  for(nChr in arrChrs) {
    
    nCount <- sum(as.integer(datConsurfRec9[datConsurfRec9$chr.x ==  nChr, sHap1]  == 1 ), na.rm=T)
    datPerHapPerChrDel <- rbind(datPerHapPerChrDel, data.frame(hap1=sHap1, chr=nChr, DelCount=nCount));
    
  }
}

for(nHap1 in 1:(length(arrHapNames)-1) ) {
  sHap1 <- arrHapNames[nHap1];
  
  
  
  for(nHap2 in (nHap1+1):length(arrHapNames) ) {
    sHap2 <- arrHapNames[nHap2];
    
    cat("\nbetween ", sHap1, " ", sHap2 ,"\n\n");
    
    nBothOK <- sum(as.integer(datConsurfRec9[ , sHap1]  == 0 & datConsurfRec9[ , sHap2] == 0))
    nHap2Del <- sum(as.integer(datConsurfRec9[ , sHap1]  == 0 & datConsurfRec9[ , sHap2] == 1))
    nHap1Del <- sum(as.integer(datConsurfRec9[ , sHap1]  == 1 & datConsurfRec9[ , sHap2] == 0))
    nBothHapDel <- sum(as.integer(datConsurfRec9[ , sHap1]  == 1 & datConsurfRec9[ , sHap2] == 1))
    datHyperGeoTests <- rbind(datHyperGeoTests, data.frame(hap1=sHap1, hap2=sHap2, BothOK=nBothOK, hap1onlydel=nHap1Del, hap2onlydel=nHap2Del, bothhapdel=nBothHapDel,
                                                           hyperP=phyper(nBothHapDel, nHap2Del + nBothHapDel, nHap1Del + nBothOK, nHap1Del + nBothHapDel, lower.tail = FALSE)));
    
    for(nChr in arrChrs) {
      cat("\nbetween ", sHap1, " ", sHap2 , " chr " , nChr, "\n\n");
      
      nBothOK <- sum(as.integer(datConsurfRec9[datConsurfRec9$chr.x ==  nChr, sHap1]  == 0 & datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap2] == 0))
      nHap2Del <- sum(as.integer(datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap1]  == 0 & datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap2] == 1))
      nHap1Del <- sum(as.integer(datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap1]  == 1 & datConsurfRec9[ datConsurfRec9$chr.x ==  nChr , sHap2] == 0))
      nBothHapDel <- sum(as.integer(datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap1]  == 1 & datConsurfRec9[datConsurfRec9$chr.x ==  nChr , sHap2] == 1))
      datHyperGeoTestsPerChr <- rbind(datHyperGeoTestsPerChr, data.frame(hap1=sHap1, hap2=sHap2, chr=nChr, BothOK=nBothOK, hap1onlydel=nHap1Del, hap2onlydel=nHap2Del, bothhapdel=nBothHapDel,
                                                                         hyperP=phyper(nBothHapDel, nHap2Del + nBothHapDel, nHap1Del + nBothOK, nHap1Del + nBothHapDel, lower.tail = FALSE)));
      
    }
  }
}



datHyperGeoTests$hap1overlappercent <- datHyperGeoTests$bothhapdel / (datHyperGeoTests$hap1onlydel + datHyperGeoTests$bothhapdel)
datHyperGeoTests$hap2overlappercent <- datHyperGeoTests$bothhapdel / (datHyperGeoTests$hap2onlydel + datHyperGeoTests$bothhapdel)
write.table(datHyperGeoTests, file="hypergeometric_test_overlap_hap_consurfdelvars.txt", col.names = T, row.names = F, quote = F, sep="\t");

datHyperGeoTestsPerChr$hap1overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap1onlydel + datHyperGeoTestsPerChr$bothhapdel)
datHyperGeoTestsPerChr$hap2overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap2onlydel + datHyperGeoTestsPerChr$bothhapdel)
write.table(datHyperGeoTestsPerChr, file="hypergeometric_test_overlap_hap_consurfdelvars_perchr.txt", col.names = T, row.names = F, quote = F, sep="\t");

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
  colnames(datac)[which(colnames(datac)=="mean")] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


datHapDelPerChrPlot <- summarySE(datHyperGeoTestsPerChrPlot, measurevar="count", groupvars=c("chr","type"))
ggplot(datHapDelPerChrPlot, aes(fill=type, y=count, x=chr)) + 
  geom_bar(position='dodge',  stat="identity") +
  geom_errorbar(aes(ymin=count-sd,
                    ymax=count+sd),
                size=0.2, width=0.7, linetype=1, position=position_dodge(.9)) +
  theme_classic() 
ggsave("consurfdel_sd_per_chr.pdf", width=12, height=2.5)


sChrFile <- "../hifiasm_haplotype_gene_variants/GOanalysis/chromosomes.tsv";


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

#
datPresence <- read.table("../hifiasm_haplotype_gene_variants/counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)

datPerHapPerChrDel <- merge(datPerHapPerChrDel, datChr, by.x="chr", by.y="Chrom", all.x = T, all.y =F);
summary(lm(datPerHapPerChrDel$DelCount ~ datPerHapPerChrDel$End))
arrGeneCountsPerChr <- table(datPresence$pgChr)[1:24]
datGeneCountsPerChr <- data.frame(chr=names(arrGeneCountsPerChr), genecount=arrGeneCountsPerChr )
datPerHapPerChrDel <- merge(datPerHapPerChrDel, datGeneCountsPerChr, by.x="chr", by.y="chr", all.x = T, all.y =F);
summary(lm(datPerHapPerChrDel$DelCount ~ datPerHapPerChrDel$genecount.Freq ))

datDelCountFoldDiffHaps <- aggregate(datPerHapPerChrDel$DelCount ,by = list(chr=datPerHapPerChrDel$chr) , FUN=function(x) {max(x)/min(x)})

arrHaplotypeIDs <- arrHapNames
datROHSegments <- NULL;
arrColors <- (as.character(wes_palette("Zissou1", length(arrHaplotypeIDs), type = "continuous")));
arrColors <- paste0(arrColors, 'aa');
nColor <- 1;
for (sF in arrHaplotypeIDs ) {
  
  datThis <- datConsurfRec9[datConsurfRec9[,sF] == 1, 10:12 ]
  colnames(datThis) <- c("chr", "start", "end");
  datThis$Group <- sF;
  datThis$Colors <- arrColors[nColor]
  nColor <- nColor + 1;
  datROHSegments <- rbind(datROHSegments , datThis);
}

datROHSegments$len <- datROHSegments$end - datROHSegments$start + 1
arrTooSmallWin <- datROHSegments$len<nMinWinPlotLen;
datROHSegments$start[arrTooSmallWin] <- datROHSegments$start[arrTooSmallWin] - (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2
datROHSegments$end[arrTooSmallWin] <- datROHSegments$end[arrTooSmallWin] + (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2

datROHSegmentsForPlot <- data.frame(Chrom=gsub("bahahascf_","", datROHSegments$chr), Start=datROHSegments$start, End=datROHSegments$end, Group=factor(x = datROHSegments$Group, levels = unique(datROHSegments$Group)), Colors=factor(datROHSegments$Colors, levels = arrColors) )

datROHSegmentsForPlot <- datROHSegmentsForPlot[as.integer(datROHSegmentsForPlot$Chrom)<=24, ]


datGeneLabels <- datConsurfRec9[, c('chr.y', 'start', 'end', 'GeneSymbol')]
datGeneLabels <- datGeneLabels[(!grepl('unknown', datGeneLabels$GeneSymbol, fixed = T)) & (!grepl('loc[0-9]*',datGeneLabels$GeneSymbol )) & (!grepl('si:',datGeneLabels$GeneSymbol )), ];
datGeneLabelsPlot <- data.frame(Chrom=gsub("bahahascf_","", datGeneLabels$chr.y), Start=datGeneLabels$start, End=datGeneLabels$end, ID=datGeneLabels$GeneSymbol );

datGeneLabelsPlot <- datGeneLabelsPlot[!duplicated(datGeneLabelsPlot),]

datROHSegmentsForPlot <- rbind(datROHSegmentsForPlot, datFakeMarkers)

pdf(file = "consurf.chr.del.pdf",width=10/1.5,height=8/1.5);
chromPlot(gaps=datChr, bands=datROHSegmentsForPlot, figCols=12, stat=datGeneLabelsPlot , 
          statCol="Value",  
          noHist = T,
          statName ="genes", statTyp="n", statSumm="none" )




dev.off();
