#setwd("D:\\raycui\\projects\\Macropodus_hongkongensis\\gatk\\slidingwin_het");
args <- commandArgs(T);
sOut <- args[2]; #"mhk.hom.mid.win.tsv"
arrChr <- paste0("bahahascf_" , 1:24);
datHet <- read.table( gzfile(args[1],'rt') , header=F);
#hist(log10(datHet$V9), breaks=1000)
nHomCutoff <- 1e-05 #if at this cutoff, count as homozygous
nHomCutoff2 <- 2e-5 ; #if at this cutoff, count as ambiguous
nMinEffSites <- 10000; #count as ambiguous if lower than this number of effective sites
nMaxAmbgStretch <- 1000000; #max ambiguous stretch before breaking block.

datHomWin <- NULL;
write.table(data.frame(chr='chr', ambstart='ambstart', homstart='homstart',  homend='homend', ambend='ambend', stringsAsFactors = F), file = sOut , row.names = F, col.names = F, sep="\t", quote = F);
for(sChr in arrChr) {
  datHetChr <- datHet[datHet$V1 == sChr, ];

  nPrevState <- -1;
  nAmbStart <- -1;
  nAmbEnd <- -1;
  nHomStart <- -1;
  nHomEnd <- -1;
  nPrevEnd <- -1;
  nPrevStart <- -1;
  nPrevMid <- -1;
  for(i in 1:nrow(datHetChr)) {
    nPosLeft <- datHetChr[i, 2];
    nPosRight <-datHetChr[i, 3];
    nPosMid <- datHetChr[i, 4];
    nEff <- datHetChr[i, 5];
    nTheta <- datHetChr[i, 9];
    nState <- 1; #het 
    if (nEff < nMinEffSites || (nTheta <=nHomCutoff2 && nTheta > nHomCutoff) ) {
      nState <- 2;
    } else if (nTheta <= nHomCutoff) {
      nState <- 0;
    }
    
    bBreakAmb <- F;

    if (nPrevState == -1) {
      nPrevState <- nState;
      nPrevStart <- 1;
      nPrevEnd <- 1;
    }
    
    if (nState ==2) {
      if (nPrevState == nState) {
        if ((nPosMid - nAmbStart+1) > nMaxAmbgStretch) {
          bBreakAmb <- T;
          if (nHomStart ==-1) {
            nAmbStart <- nPosMid;
          }
        }
      } else if (nPrevState == 1) {
        nAmbStart <- nPosMid;
      } else if (nPrevState == 0) {
        nHomEnd <- nPrevEnd;
        nAmbEnd <- nPosMid;
      }
      
    }
    
    if (nState == 0) {
      if (nPrevState == nState) {
        nHomEnd <- nPosMid;
      } else if (nPrevState == 2) {
        if (nAmbStart < nHomStart) {
          nHomEnd <- nPosMid;
        } else {
          nHomStart <- nPosMid;
          nHomEnd <- nPosMid;
        }
      } else if (nPrevState == 1) {
       # cat("0 to 1:", nPosLeft, "\n");
        nAmbStart <- nPosMid;
        nHomStart <- nPosMid;
      }
    }
    
    if (nState == 1) {
      if (nPrevState == nState) {
        
      } else if (nPrevState == 2) {
        nAmbEnd <- nPrevEnd;
        bBreakAmb <- T;
      } else if (nPrevState == 0) {
        nAmbEnd <- nPrevEnd;
        nHomEnd <- nPrevEnd;
        bBreakAmb <- T;
      }
    }
    
    if (bBreakAmb) {
      if (nHomStart > 0 && nAmbEnd >= nHomEnd) {
        #datHomWin <- rbind(datHomWin, data.frame(chr=sChr, ambstart=nAmbStart, homstart=nHomStart,  homend=nHomEnd, ambend=nAmbEnd, stringsAsFactors = F));
        if (nAmbStart == -1) {
          nAmbStart <- nHomStart;
        }
        write.table(data.frame(chr=sChr, ambstart=nAmbStart, homstart=nHomStart,  homend=nHomEnd, ambend=nAmbEnd, stringsAsFactors = F), file = sOut , append = T, row.names = F, col.names = F, sep="\t", quote = F);
        
      }
      if (nAmbStart != nPosMid) {
      nAmbStart <- -1;
      }
      nAmbEnd <- -1;
      nHomStart <- -1;
      nHomEnd <- -1;
      
    }
    nPrevState <- nState;
    nPrevEnd <- nPosRight;
    nPrevStart <- nPosLeft;
    nPrevMid <- nPosMid;
  }
}
