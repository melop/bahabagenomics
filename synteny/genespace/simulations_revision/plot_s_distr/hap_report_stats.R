#setwd("/data/projects/rcui/bahaha_assembly/synteny/genespace/haplotype_gene_variants/simulations/");

nChr <- 24
datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)
datChrLens <- read.table("../../../../annotations/repeatmask/female.hap0/adjust_chr_order/hap.F0.fa.fai", header=F, stringsAsFactors = F, sep="\t")
arrChrLens <- datChrLens[1:nChr, 2]
arrChrLens <- arrChrLens / sum(arrChrLens);

datP <- datPresence[,5:ncol(datPresence)]
datH <- datHighImpact[,5:ncol(datHighImpact)]

nSamples <- ncol(datP)
datRet <- datH

for(nIdx in 1:ncol(datP) ) {
  datRet[, nIdx] <- as.integer(datP[ , nIdx] >=0.5 & datH[, nIdx]==0)
}


datOnlyComplete <- cbind(datPresence, datH);

datOnlyComplete <- datOnlyComplete[ rowSums( datOnlyComplete[,5:(5+nSamples-1)])==4 , ]

datOnlyCompHighImpact <- datOnlyComplete[, (5+nSamples):ncol(datOnlyComplete)]

datOnlyCompHighImpact[datOnlyCompHighImpact>0] <- 1;
datOnlyCompHighImpact <- cbind(datOnlyComplete[,1:4], datOnlyCompHighImpact )

datOnlyCompHighImpact <- datOnlyCompHighImpact[datOnlyCompHighImpact$pgChr %in% 1:nChr, ];

#get haplotype stats
datHaps <- datOnlyCompHighImpact[, 5:8]
arrHomDel <- c();
arrHet <- c();
arrDel <- c();
for(n1 in 1:(ncol(datHaps)-1)) {
  for(n2 in 1:(ncol(datHaps))) {
    if (n1>=n2) {next;}
    datTestHaps <- datHaps[, c(n1,n2)];
    arrStates <- rowSums(datTestHaps)
    arrHomDel <- c(arrHomDel, sum(as.integer(arrStates==2)));
    arrHet <- c(arrHet, sum(as.integer(arrStates==1)));
    arrDel <- c(arrDel, colSums(datTestHaps));
  }
}

#create parental populations:
nMeanHomDel <- mean(arrHomDel);
nSDHomDel <- sd(arrHomDel);
nMeanDel <- mean(arrDel);
nSDDel <- sd(arrDel);

datInfo <- data.frame(MeanHomDel = nMeanHomDel, SDHomDel=nSDHomDel, MeanDel = nMeanDel, SDDel = nSDDel );
write.table(datInfo, file="hap_del_stats.txt", row.names=F, col.names=T, sep="\t", quote=F);

fnGenerateHap <- function() {
  arrRet <- rep(0, nrow(datHaps));
  nDel <- round( rnorm(1, nMeanDel, nSDDel) );
  
  for(nHap in 1:ncol(datHaps)) {
    arrRef <- datHaps[ , nHap];
    nHomDel <- round( rnorm(1, nMeanHomDel, nSDHomDel) );
    arrDelPos <- which(arrRef == 1);
    arrSetDelPos <- sample(arrDelPos, nHomDel, replace = F) 
    arrRet[arrSetDelPos] <- 1;
  }
  
  nLeft <- nDel - sum(arrRet);
  if (nLeft <=0) {
    nLeft <- 1;
  }
  
  arrOccupiedPos <- which(arrRet==1);
  arrPickPos <- 1:length(arrRet);
  arrPickPos <- sample(arrPickPos[!(arrPickPos %in% arrOccupiedPos) ], nLeft, replace = F );
  arrRet[arrPickPos] <- 1;
  return(arrRet);
}

