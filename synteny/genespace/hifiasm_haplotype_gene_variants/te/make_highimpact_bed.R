setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/te");
library(chromPlot)
library(RColorBrewer)
library(wesanderson)

arrChrs <- 1:24;
nMinWinPlotLen <- 3e5;
sRefSp <- "BahahaTaipingensis"

datRNACoord <- read.table("../GOanalysis/mRNA.coord.txt", header=T, sep="\t", quote="")
datOrthID2RNAID <- read.table("../../hifiasm_bahaha_haplotypes/synorthos.txt", header=T, sep="\t")
datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)

arrSampleIDs <- colnames(datHighImpact)[5:ncol(datHighImpact)];
datP <- datPresence[,5:ncol(datPresence)]
datH <- datHighImpact[,5:ncol(datHighImpact)]

nSamples <- ncol(datP)
arrSamples <- colnames(datP)
datRet <- datH

for(nIdx in 1:ncol(datP) ) {
  datRet[, nIdx] <- as.integer(datP[ , nIdx] < 0.5) + datH[, nIdx]
}


datOnlyCompHighImpact <- cbind(datPresence[,1:3], datRet);
datOnlyCompHighImpact$UniqID <- paste(datOnlyCompHighImpact$pgChr, datOnlyCompHighImpact$pgOrd, datOnlyCompHighImpact$pgID, sep="_")
datOrthID2RNAID$UniqID <- paste(datOrthID2RNAID$pgChr, datOrthID2RNAID$pgOrd, datOrthID2RNAID$pgID, sep="_")
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datOrthID2RNAID[, c("UniqID",  sRefSp)], by="UniqID", all.x = T, all.y=F);

arrNonRepresentativeInDupArray <- datOnlyCompHighImpact$BahahaTaipingensis[! (datOnlyCompHighImpact$BahahaTaipingensis %in% datRNACoord$mRNAID) ]

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datRNACoord, by.x = sRefSp, by.y="mRNAID");
datOnlyCompHighImpact$all_haps_mean_del <- rowSums(datOnlyCompHighImpact[, arrSampleIDs])/length(arrSampleIDs)

for(sSample in c(arrSampleIDs, 'all_haps_mean_del') ) {
  datBed <- datOnlyCompHighImpact[ datOnlyCompHighImpact[, sSample] >0, c('chr', 'start', 'end', sSample)];
  datBed$mid <- as.integer( (datBed$end + datBed$start)/2)
  datBed$start <- datBed$mid-1;
  datBed$end <- datBed$mid+ as.integer(datBed[,sSample]);
  write.table(datBed[, c('chr', 'start', 'end')], file=paste0("highimpact.",sSample,".bed"), col.names = F, row.names = F, quote=F, sep="\t");
}

#system("bash getrhooverlap.sh")
