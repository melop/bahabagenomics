#setwd("/data/projects/rcui/bahaha_assembly/synteny/genespace/haplotype_gene_variants/simulations/");
args <- commandArgs(trailingOnly=T);

nSeed <- as.integer(args[1]); #338823232;
cat("seed is ", nSeed, "\n");
set.seed(nSeed)  
nChr <- as.integer(args[2]); #24;
nTotalGen <- as.integer(args[3]); #100;
nSeedingPopSize <-  as.integer(args[4]); # 10;#starting population size
nPopSize <- as.integer(args[5]); # 10;
h <- as.numeric(args[6]); #  0.5; 
s <- as.numeric(args[7]); # 0.002;#0.001; #selection coefficient per locus model following w(i) = max( 1 - h * s * i1 - s*i2, 0); i1 is number het deleterious sites, i2 is number of hom del sites
nRecBreakpointPerChr <-as.numeric(args[8]); #  20; #recombination breakpoints per chromosome per gamete
nKeepLoadPercentile <- as.numeric(args[9]); # per generation only keep this percentile of individuals that have lowest load per genome

nWildNe <- 0.00131972/(4*7.75273E-9); #theta from ~/bahaha_assembly/wgs/gatk/slidingwin_het_100kb_10kb/, mu from ~/bahaha_assembly/wgs/psmc/called_separately/plotting.R 


cat("nChr = ", nChr, "\nnTotalGen = ", nTotalGen, "\nnSeedingPopSize=",nSeedingPopSize, "\nnPopSize=",nPopSize,"\nh = ", h, "\ns = ", s, "\nnRecBreakpointPerChr = ",nRecBreakpointPerChr,"\nnKeepLoadPercentile=",nKeepLoadPercentile,"\n");
nMaxTryPerTime <- nPopSize;

# # Calculate the number of cores
# no_cores <- 24;
# 
# # Initiate cluster
# cl <- makeCluster(no_cores)


datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datHighImpact <- read.table("../counts_deleterious_ret.HIGH.txt", header = T, stringsAsFactors = F)
datChrLens <- read.table("../../../../annotations/repeatmask/female.hap0/adjust_chr_order/hap.F0.fa.fai", header=F, stringsAsFactors = F, sep="\t")
arrChrLens <- datChrLens[1:nChr, 2]
arrChrLens <- arrChrLens / sum(arrChrLens);
arrRecBpPerChr <- nRecBreakpointPerChr * nChr * arrChrLens

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

#compute mutation probability per haplotype per meiosis
fnWattersonThetaDenominator <- function(n) {
  nSum <- 0;
  for(i in 2:n) {
    nSum <- nSum + 1/(i-1);
  }
  return(nSum)
}

arrStates <- rowSums(datHaps);
thetaWatterson <- sum(as.integer(arrStates>0 & arrStates<ncol(datHaps) ))/fnWattersonThetaDenominator(ncol(datHaps))
U <- thetaWatterson / (4*nWildNe);

cat("New deleterious mutations per haplotype per generation : ", U, "\n");

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

lsParentalPop <- lapply(1:nSeedingPopSize, function(x) {
    datInd <- datOnlyCompHighImpact[, c(1,5,6)];
    colnames(datInd)[2:3] <- c('hap1','hap2');
    datInd[, 2] <- fnGenerateHap();
    datInd[, 3] <- fnGenerateHap();
    return(datInd);
  }
  )

fnFitness <- function(nHet, nHomDel) {
  return(max(1- h * s* nHet - s * nHomDel, 0));
}

# sum(datOnlyCompHighImpact$BahahaTaipingensisF0)/nrow(datOnlyCompHighImpact);
# sum(datOnlyCompHighImpact$BahahaTaipingensisF1)/nrow(datOnlyCompHighImpact);
# sum(datOnlyCompHighImpact$BahahaTaipingensisM0)/nrow(datOnlyCompHighImpact);
# sum(datOnlyCompHighImpact$BahahaTaipingensisM1)/nrow(datOnlyCompHighImpact);

#fnFitness(sum(datOnlyCompHighImpact$BahahaTaipingensisF0));
lsPop <- lsParentalPop;
#View(lsPop[[2]])
nKilledInd <- 0;
nTotalTries <-0;
fnBreed <- function(lsPop) {
  nInd <- length(lsPop);
  lsNewPop <- list();
  if (nInd <= 1) {
    cat("Not enough individuals.");
    return(NA)
  }
  
  nKilledInd <<- 0;
  nTotalTries <<-0;
  arrFitnessBeforeSel <<- c();
  nTryTimes <- nPopSize;
  while (T) {
    
    lsThisTry <- lapply(1:nTryTimes, function (nRep) {
      
      nParent1 <- sample.int(nInd, 1)
      nParent2 <- -1;
      while (nParent1==nParent2 || nParent2==-1 || ((nParent1 %% 2) == (nParent2 %% 2) ) ) { #must be different sexes to breed
        nParent2 <- sample.int(nInd, 1);
      }
      
      nTotalRows <- nrow(lsPop[[nParent1]]);
      datOffspring <- data.frame(pgChr=integer(nTotalRows), hap1=integer(nTotalRows), hap2=integer(nTotalRows) ); 
      datChr1 <- fnGetGamete(lsPop[[nParent1]]);
      datChr2 <- fnGetGamete(lsPop[[nParent2]]);
      datOffspring[,1] <- datChr1[,1];
      datOffspring[,2] <- datChr1[,2];
      datOffspring[,3] <- datChr2[,2];
      
      #now compute fitness
      arrSiteTypes <- rowSums(datOffspring[,2:3]);
      nHet <- sum(as.integer(arrSiteTypes == 1));
      nHomDel <- sum(as.integer(arrSiteTypes == 2));
      w <- fnFitness(nHet, nHomDel);
      arrFitnessBeforeSel <<- c(arrFitnessBeforeSel, w);
      nRand <- runif(n = 1, min = 0, max = 1)
      if (nRand <= w) { #keep individual
        return(datOffspring);
      } else {
        #nKilledInd <<- nKilledInd + 1;
        return(NA);
      }
    })
    
    nTotalTries <<- nTotalTries + nTryTimes;
    nKilled <- sum(as.integer(is.na(lsThisTry)));
    nKilledInd <<- nKilledInd + nKilled;
    
    nKilledPerc <- nKilledInd / nTotalTries;
    
    lsThisTry <- lsThisTry[!is.na(lsThisTry)]; #remove dead individuals
    arrLoadSum <- unlist( lapply(lsThisTry, function(x) {sum(x)} ) );
    nCutoff <- quantile(arrLoadSum, nKeepLoadPercentile);
    arrKeep <- which(arrLoadSum <= nCutoff);
    lsThisTry <- lsThisTry[arrKeep];
    
    nStillNeeded <- nPopSize - length(lsNewPop);
    nCopyInds <- 1;
    if (nStillNeeded > length(lsThisTry)) {
      nCopyInds <- length(lsThisTry);
      nAdjustFactor <- (1 - nKilledPerc)/nKeepLoadPercentile;
      if (nAdjustFactor<0) {
        cat("Warning: selection coefficient may be too high to produce viable offspring.\n");
        nAdjustFactor <- 0.5;
      }
      nTryTimes <- ceiling( (nStillNeeded - nCopyInds )/nAdjustFactor);
      if (nTryTimes>nMaxTryPerTime) {
        nTryTimes<-nMaxTryPerTime;
      }
    } else {
      nCopyInds <- nStillNeeded;
      nTryTimes <- 0;
    }
    
    if (nCopyInds>=1) {
      lsNewPop <- c(lsNewPop, lsThisTry[1:nCopyInds]);
    }
    if (nTryTimes<=0) {
      break;
    }
    
  }
  
  return(lsNewPop);
  
}

fnPopStats <- function(lsPop,nGen) {
    arrHet <- c();
    arrHomDel <- c();
    arrFitness <- c();
    arrHapDelCount <- c();
    for(nInd in 1:length(lsPop)) {
      datInd <- lsPop[[nInd]];
      arrSiteTypes <- rowSums(datInd[,2:3]);
      nHet <- sum(as.integer(arrSiteTypes == 1));
      nHomDel <- sum(as.integer(arrSiteTypes == 2));
      w <- fnFitness(nHet, nHomDel);
      arrHapDelCount <- c(arrHapDelCount , colSums(datInd[, 2:3]));
      
      arrHet <- c(arrHet, nHet/length(arrSiteTypes));
      arrHomDel <- c(arrHomDel, nHomDel/length(arrSiteTypes));
      arrFitness <- c(arrFitness, w);
      
    }
    
    return(data.frame(gen=nGen, 
                      het_mean=mean(arrHet),
                      het_sd=sd(arrHet),
                      homdel_mean=mean(arrHomDel),
                      homdel_sd=sd(arrHomDel),
                      fitness_beforesel_mean=mean(arrFitnessBeforeSel),
                      fitness_beforesel_sd=sd(arrFitnessBeforeSel),
                      fitness_postsel_mean=mean(arrFitness),
                      fitness_postsel_sd=sd(arrFitness),
                      delperhap_mean=mean(arrHapDelCount),
                      delperhap_sd=sd(arrHapDelCount),
                      selection_killedind=nKilledInd / nTotalTries
                      ));
}


fnGetGamete <- function(datInd) {
  datGamete <- datInd;
  colnames(datGamete)[2:3] <- c("hap1", "hap2");
  datGamete$new <- 0;
  for (nThisChr in 1:nChr) {
    nBreakPoints <- rpois(1,arrRecBpPerChr[nChr]);
    nFirstChr <- sample.int(2,1);
    nSecondChr <- 0;
    if (nFirstChr == 1) {
      nSecondChr <- 2;
    } else {
      nSecondChr <- 1;
    }
    
    
    nFirstChr <- nFirstChr + 1;
    nSecondChr <- nSecondChr + 1;
    datChr <- datGamete[datGamete$pgChr==nThisChr, ];
    
    if (nBreakPoints <=0) {
      datGamete[datGamete$pgChr==nThisChr, 'new' ] <- datChr[ , nFirstChr ] ;
      next;
    }
    

    arrBps <- sort(c( sample.int(nrow(datChr)-1,nBreakPoints ),nrow(datChr)) );
    arrStartPts <- c(1, arrBps[1:(length(arrBps)-1)]+1 )
    bSwitch <- T;
    for(nBlock in 1:length(arrStartPts)) {
      nStartRow <- arrStartPts[nBlock];
      nEndRow <- arrBps[nBlock];
      nThisHap <- nFirstChr;
      if (!bSwitch) {
        nThisHap <- nSecondChr;
      }
      bSwitch <- (!bSwitch);
      
      datChr[nStartRow:nEndRow, 'new'] <- datChr[nStartRow:nEndRow, nThisHap];
    }
    
    nNewDelMut <- rpois(1, U);
    if (nNewDelMut > 0) {
      cat("New del mut:", nNewDelMut, "\n");
      arrMutPos <- sample(which(datChr[, 'new']!=1), nNewDelMut);
      datChr[arrMutPos , 'new'] <- 1;
    }
    datGamete[datGamete$pgChr==nThisChr, 'new' ] <- datChr[, 'new'];
  }
  
  return(datGamete[, c('pgChr', 'new')]);
}

arrFitnessBeforeSel <-1;

#clusterExport(cl, c("fnGetGamete", "nChr", "nRecBreakpointPerChr", "fnFitness", 'h', 's','arrFitnessBeforeSel' ) )

datSimStats <- fnPopStats(lsPop,0);
for(nGen in 1:nTotalGen) {
  cat(nGen,"\n")
  lsPop <- fnBreed(lsPop);
  datSimStats <- rbind(datSimStats, fnPopStats(lsPop,nGen));
  
}

datPlot <- datSimStats[2:nrow(datSimStats), ]
sOutStem <- paste("simret","gen",nTotalGen,"seedpopsize",nSeedingPopSize, "popsize",nPopSize, "h",h, "s", s, 'reccount',nRecBreakpointPerChr, 'nKeepLoadPercentile', nKeepLoadPercentile, 'seed', nSeed, sep="_")

write.table(datSimStats, file = paste0(sOutStem, ".txt"), quote = F, row.names = F, col.names = T, sep="\t")

#plot(datPlot$gen, datPlot$selection_killedind)
