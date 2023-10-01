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

#summary(lm(datOnlyCompHighImpact$BahahaTaipingensisF0 ~ datOnlyCompHighImpact$BahahaTaipingensisF1))

#arrHapNames <- colnames(datOnlyCompHighImpact)[6:7] #female
arrHapNames <- colnames(datOnlyCompHighImpact)[8:9] # male

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
#write.table(datHyperGeoTests, file="hypergeometric_test_overlap_hap_delvars_female.txt", col.names = T, row.names = F, quote = F, sep="\t");
write.table(datHyperGeoTests, file="hypergeometric_test_overlap_hap_delvars_male.txt", col.names = T, row.names = F, quote = F, sep="\t");

datHyperGeoTestsPerChr$hap1overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap1onlydel + datHyperGeoTestsPerChr$bothhapdel)
datHyperGeoTestsPerChr$hap2overlappercent <- datHyperGeoTestsPerChr$bothhapdel / (datHyperGeoTestsPerChr$hap2onlydel + datHyperGeoTestsPerChr$bothhapdel)
#write.table(datHyperGeoTestsPerChr, file="hypergeometric_test_overlap_hap_delvars_perchr_female.txt", col.names = T, row.names = F, quote = F, sep="\t");
write.table(datHyperGeoTestsPerChr, file="hypergeometric_test_overlap_hap_delvars_perchr_male.txt", col.names = T, row.names = F, quote = F, sep="\t");

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




ggplot(datHyperGeoTestsPerChrPlot, aes(fill=type, y=count, x=chr)) + 
  geom_bar(position='dodge',  stat="identity") +
   theme_classic() 
#ggsave("del_sd_per_chr_female.pdf", width=12, height=2.5)
ggsave("del_sd_per_chr_male.pdf", width=12, height=2.5)

datFemale <- read.table("hypergeometric_test_overlap_hap_delvars_perchr_female.txt", sep="\t", header=T);
datMale <- read.table("hypergeometric_test_overlap_hap_delvars_perchr_male.txt", sep="\t", header=T);
datFemale$meanshared <- (datFemale$hap1overlappercent + datFemale$hap2overlappercent)/2
datMale$meanshared <- (datMale$hap1overlappercent + datMale$hap2overlappercent)/2
wilcox.test(datFemale$meanshared, datMale$meanshared, paired = T)
mean(datFemale$meanshared);
mean(datMale$meanshared);

wilcox.test(datFemale$bothhapdel, datMale$bothhapdel, paired = T)
mean(datFemale$bothhapdel - datMale$bothhapdel);





datFemale$delsites <- datFemale$hap1onlydel + datFemale$hap2onlydel + datFemale$bothhapdel
datMale$delsites <- datMale$hap1onlydel + datMale$hap2onlydel + datMale$bothhapdel
wilcox.test(datFemale$delsites,  datMale$delsites, paired = T)


datDiffPlot <- data.frame(chr=as.integer(), diff=as.integer(), type=character());
datAdd <- data.frame(chr=datFemale$chr, diff=(datFemale$hap1onlydel + datFemale$hap2onlydel + 2*datFemale$bothhapdel)/2 - (datMale$hap1onlydel + datMale$hap2onlydel + 2*datMale$bothhapdel)/2
                     , type="1del_per_hap_diff");
datDiffPlot <- rbind(datDiffPlot, datAdd);
datAdd <- data.frame(chr=datFemale$chr, diff=(datFemale$hap1onlydel + datFemale$hap2onlydel)/2  - (datMale$hap1onlydel + datMale$hap2onlydel)/2 
                     , type="2private_to_hap_diff");
datDiffPlot <- rbind(datDiffPlot, datAdd);
datAdd <- data.frame(chr=datFemale$chr, diff=datFemale$bothhapdel - datMale$bothhapdel
                     , type="3bothhapdel_diff");
datDiffPlot <- rbind(datDiffPlot, datAdd);


ggplot(datDiffPlot, aes(fill=type, y=diff, x=chr)) + 
  geom_bar(position='dodge',  stat="identity") +
  theme_classic() 
ggsave("del_sd_per_chr_diff_female-male.pdf", width=12, height=2.5)

