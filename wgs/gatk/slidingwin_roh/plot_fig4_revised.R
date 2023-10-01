setwd("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_roh");
library(chromPlot)
library(RColorBrewer)
library(wesanderson)
library(svglite)

#pdf(file = "fig4.roh.BTP_Q0.revised.pdf",width=10,height=8);
#pdf(file = "fig4.roh.parentalpop.revised.pdf",width=10,height=8);
svglite(file = "fig4.roh.parentalpop.revised.svg",width=10,height=8);
#svglite(file = "fig4.roh.BTP_Q0.revised.svg",width=10,height=8);
#svglite(file = "fig4.roh.chrbackground.revised.svg",width=10,height=8);

nMinWinPlotLen <- 3e5;

sChrFile <- "/data/projects/rcui/bahaha_assembly/wgs/sex_chrm_id_q0q1/chromosomes.tsv";

#roh histogram
datROH <- read.table("roh.tsv", sep="\t", header=F)

nNGS_Longread_factor <- datROH[ datROH$V1 == 'BTP_Q1 ', 9] /datROH[ datROH$V1 == 'BTP_M06 ', 9]

datROH$calibratedFROH <- datROH$V9 * nNGS_Longread_factor;
datROH$calibratedFROH[datROH$V1 == 'BTP_Q0 '] <- datROH$V9[datROH$V1 == 'BTP_Q0 '];#third gen data
datROH$calibratedFROH[datROH$V1 == 'BTP_Q1 '] <- datROH$V9[datROH$V1 == 'BTP_Q1 '];#third gen data

datROH <- datROH[datROH$V1!='BTP_Q1 ', ]
hist(datROH$calibratedFROH, breaks = seq(0,0.12, 0.01), col = "#3B9AB2" )

write.table(datROH, file="roh.calibrated.tsv", sep="\t", col.names = F, row.names = F, quote=F)

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
datFakeMarkers$Colors <- '#D3D3D3'
  
arrHomFiles <- Sys.glob("*.hom.mid.win.tsv");
#remove M06 since it is same individual as Q1
arrHomFiles <- arrHomFiles[!grepl("BTP_Q1", arrHomFiles)]
#arrHomFiles <- arrHomFiles[grepl("BTP_Q0", arrHomFiles)]
arrHomFiles <- arrHomFiles[!grepl("BTP_Q0", arrHomFiles)]

length(arrHomFiles)
datROHSegments <- NULL;
arrColors <- as.character(wes_palette("Zissou1", length(arrHomFiles), type = "continuous"))
#arrColors <- paste0(arrColors, '99');
nColor <- 1;
for (sF in arrHomFiles) {
  datThis <- read.table(sF, header=T, sep="\t");
  datThis$Group <- gsub('BTP_','', gsub('.hom.mid.win.tsv','', sF, fixed=T), fixed=T);
  datThis$Colors <- arrColors[nColor]
  nColor <- nColor + 1;
  datROHSegments <- rbind(datROHSegments , datThis);
}



#recombination rate:
datRec <- read.table("/data/projects/rcui/bahaha_assembly/ismc/20ind/20ind/rho.excludeROH.txt", header=F, sep="\t") 
datRec$mid <- (datRec$V2 + datRec$V3)/2
datROHSegments$recrate <- 0;
arrTargetRecRate <- NULL;
for (nRow in 1:nrow(datROHSegments) ) {
  arrRecs <-datRec[datRec$V1 ==  datROHSegments[nRow, 'chr'] & datRec$mid >= datROHSegments[nRow, 'homstart'] & datRec$mid <=  datROHSegments[nRow, 'homend'], 'V4'];
  arrTargetRecRate <- c(arrTargetRecRate, arrRecs);
}

datRecRatePlot <- data.frame(recrate = arrTargetRecRate, regiontype='roh' , stringsAsFactors = F);
datRecRatePlot <- rbind(datRecRatePlot, data.frame(recrate = datRec$V4, regiontype='control' , stringsAsFactors = F));
boxplot(log10(datRecRatePlot$recrate)~datRecRatePlot$regiontype, col = "#7EB8BC", pch="16")
wilcox.test(arrTargetRecRate, datRec$V4)
mean(arrTargetRecRate);
mean(datRec$V4);
median(arrTargetRecRate);
median(datRec$V4);

#

datROHSegments$len <- datROHSegments$homend - datROHSegments$homstart + 1
hist(log10(datROHSegments$len), xlim=c(4,7), breaks=50 , col="#E7C21C")
arrTooSmallWin <- datROHSegments$len<nMinWinPlotLen;
datROHSegments$homstart[arrTooSmallWin] <- datROHSegments$homstart[arrTooSmallWin] - (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2
datROHSegments$homend[arrTooSmallWin] <- datROHSegments$homend[arrTooSmallWin] + (datROHSegments$len[arrTooSmallWin]-nMinWinPlotLen)/2

datROHSegmentsForPlot <- data.frame(Chrom=gsub("bahahascf_","", datROHSegments$chr), Start=datROHSegments$homstart, End=datROHSegments$homend, Group=datROHSegments$Group, Colors=datROHSegments$Colors )

#add end markers to make sure chromosome length is correct:

datROHSegmentsForPlot <- rbind( datROHSegmentsForPlot, datFakeMarkers)

chromPlot( bands=datROHSegmentsForPlot, gaps=datChr, figCols=12 )
#chromPlot( gaps=datChr, figCols=12 ) #empty plot



dev.off();
