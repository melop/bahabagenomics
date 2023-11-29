#setwd("/data/projects/rcui/bahaha_assembly/synteny/genespace/haplotype_gene_variants/simulations/");

# # Calculate the number of cores
# no_cores <- 24;
# 
# # Initiate cluster
# cl <- makeCluster(no_cores)

nChr <- 24
datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)

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

write.table(datHaps, file="test.txt", col.names=T, row.names=F, sep="\t", quote=F);
