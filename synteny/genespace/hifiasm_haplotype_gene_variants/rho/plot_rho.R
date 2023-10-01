setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/rho");
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
  datRet[, nIdx] <- as.integer(datP[ , nIdx] >=0.5 & datH[, nIdx]==0)
}


datOnlyCompHighImpact <- cbind(datPresence[,1:3], datRet);
datOnlyCompHighImpact$UniqID <- paste(datOnlyCompHighImpact$pgChr, datOnlyCompHighImpact$pgOrd, datOnlyCompHighImpact$pgID, sep="_")
datOrthID2RNAID$UniqID <- paste(datOrthID2RNAID$pgChr, datOrthID2RNAID$pgOrd, datOrthID2RNAID$pgID, sep="_")
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datOrthID2RNAID[, c("UniqID",  sRefSp)], by="UniqID", all.x = T, all.y=F);

arrNonRepresentativeInDupArray <- datOnlyCompHighImpact$BahahaTaipingensis[! (datOnlyCompHighImpact$BahahaTaipingensis %in% datRNACoord$mRNAID) ]

datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datRNACoord, by.x = sRefSp, by.y="mRNAID");

write.table(datOnlyCompHighImpact[, c('chr', 'start', 'end', 'UniqID')], file="_tmp.bed", col.names = F, row.names = F, quote=F, sep="\t");

system("bash getrhooverlap.sh")

datOnlyCompHighImpactRho <- read.table("SV.rho.bed", header=F, sep="\t");
datOnlyCompHighImpactRho <- aggregate(datOnlyCompHighImpactRho[,8], by=list(UniqID=datOnlyCompHighImpactRho[,4]), FUN=mean)
colnames(datOnlyCompHighImpactRho)[2] <- "rho"
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datOnlyCompHighImpactRho, by="UniqID", all.x=T, all.y=F)
datOnlyCompHighImpact$all_haps_pooled_freq <- rowSums(datOnlyCompHighImpact[, arrSampleIDs])/length(arrSampleIDs)

datBoxPlot <- NULL;
for(sSample in c(arrSampleIDs, 'all_haps_pooled_freq') ) {
  
  arrGoodGeneRho <- datOnlyCompHighImpact[datOnlyCompHighImpact[,sSample]==1, 'rho' ]
  arrBadGeneRho <- datOnlyCompHighImpact[datOnlyCompHighImpact[,sSample]<1, 'rho' ]
  
  arrGoodGeneRho <- arrGoodGeneRho[!is.na(arrGoodGeneRho)];
  arrBadGeneRho <- arrBadGeneRho[!is.na(arrBadGeneRho)];
  cat("=============== ", sSample," =================\n");
  cat("Good genes count :",length(arrGoodGeneRho),"\n");
  cat("Bad genes  count :",length(arrBadGeneRho),"\n");
  cat("Good genes mean rho :",mean(arrGoodGeneRho),"\n");
  cat("Bad genes mean rho :",mean(arrBadGeneRho),"\n");
  cat("Good genes median rho :",median(arrGoodGeneRho),"\n");
  cat("Bad genes median rho :",median(arrBadGeneRho),"\n");
  print(wilcox.test(arrGoodGeneRho, arrBadGeneRho))
  print(t.test(arrGoodGeneRho, arrBadGeneRho))
  
 
  datBoxPlot <- rbind(datBoxPlot, data.frame(sample = sSample , rho = datOnlyCompHighImpact[, 'rho' ] , isgoodgene=(datOnlyCompHighImpact[, sSample ]==1)));
  
}

boxplot(log10(datBoxPlot$rho) ~ datBoxPlot$isgoodgene * datBoxPlot$sample );

datBoxPlot <- datBoxPlot[datBoxPlot$sample!="all_haps_pooled_freq",];
ggplot(datBoxPlot, aes(x=sample, y=log10(rho), fill=isgoodgene)) + 
  geom_boxplot()+
  theme_minimal() + 
  scale_fill_manual('High Impact Variant Type', values=rev(as.character(wes_palette("Zissou1", 2, type = "continuous"))) )

ggsave("impacttype_perhap_vs_recrate.pdf", width=6, height=6)
summary(lm(log10(rho) ~ isgoodgene + sample, data=datBoxPlot));
