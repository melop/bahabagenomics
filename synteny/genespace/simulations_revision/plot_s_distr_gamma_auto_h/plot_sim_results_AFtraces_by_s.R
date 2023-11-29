setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/simulations_revision/plot_s_distr_gamma_auto_h/");
library("wesanderson")

xlim <- c(0,200);


arrParamNameMap <- c("nchr"      ,                "gen"   ,              "N0" ,          "Ne"   ,               "h"  ,                     
                     "h_sd"      ,               "s"          ,              "s_sd"         ,            "r"   ,             "k",
                     "MeanHomDel"  ,             "SDHomDel"   ,              "MeanDel"      ,            "SDDel"     ,               "M"   ,  
                     "Reps"        ,             "OutStem"); 
names(arrParamNameMap) <- c("Chr"      ,                "TotalGen"   ,              "SeedingPopSize" ,          "PopSize"   ,               "h"  ,                     
                            "h_sd"      ,               "s"          ,              "s_sd"         ,            "RecPerChr"   ,             "TopPercentileForBreeding",
                            "MeanHomDel"  ,             "SDHomDel"   ,              "MeanDel"      ,            "SDDel"     ,               "GenotypedMarkerPerc"   ,  
                            "Reps"        ,             "OutStem");


arrYLimMax <- c(200, 1, 1, 0.15);

# ylim <- c(0,200);
# sPlotVar <- "delperhap_mean";
# sPlotVarName <- "Deleterious mutations per haplotype";
# 

# ylim <- c(0,1);
# sPlotVar <- "selection_killedind";
# sPlotVarName <- "Proportion of individuals killed due to viability selection";

# ylim <- c(0,1);
# sPlotVar <- "fitness_beforesel_mean";
# sPlotVarName <- "Mean absolute fitness";

# ylim <- c(0,0.06);
# sPlotVar <- "fitness_beforesel_sd";
# sPlotVarName <- "stddev absolute fitness";


datSimDefs <- NULL
#arrDefs <- Sys.glob("simdef.txt")
arrDefs <- c("simdef2.txt")

for(sDef in arrDefs) {
  datSimDefs <- rbind(datSimDefs, read.table(sDef, header=T, sep="\t") )
}

datSimDefs$OutStem <- paste("./out/simret","gen",datSimDefs$TotalGen,
                  "seedpopsize",datSimDefs$SeedingPopSize, 
                  "popsize",datSimDefs$PopSize, 
                  "h",datSimDefs$h, "s", datSimDefs$s, "h_sd", datSimDefs$h_sd, "s_sd", datSimDefs$s_sd,
                  'reccount',datSimDefs$RecPerChr, 
                  'nKeepLoadPercentile', datSimDefs$TopPercentileForBreeding, 
                  'nMeanHomDel',  as.integer(datSimDefs$MeanHomDel), 
                  'nSDHomDel', as.integer(datSimDefs$SDHomDel), 
                  'nMeanDel', as.integer(datSimDefs$MeanDel), 
                  'nSDDel',as.integer(datSimDefs$SDDel), 
                  'GenotypePercMarkers',datSimDefs$GenotypedMarkerPerc,
                  'seed', '*perlocusfitness.txt', sep="_")

datSimDefs$OutStem2 <- paste("./out/simret","gen",datSimDefs$TotalGen,
                            "seedpopsize",datSimDefs$SeedingPopSize, 
                            "popsize",datSimDefs$PopSize, 
                            "h",datSimDefs$h, "s", datSimDefs$s, "h_sd", datSimDefs$h_sd, "s_sd", datSimDefs$s_sd,
                            'reccount',datSimDefs$RecPerChr, 
                            'nKeepLoadPercentile', datSimDefs$TopPercentileForBreeding, 
                            'nMeanHomDel',  as.integer(datSimDefs$MeanHomDel), 
                            'nSDHomDel', as.integer(datSimDefs$SDHomDel), 
                            'nMeanDel', as.integer(datSimDefs$MeanDel), 
                            'nSDDel',as.integer(datSimDefs$SDDel), 
                            'GenotypePercMarkers',datSimDefs$GenotypedMarkerPerc,
                            'seed', '*AF.dump.tsv', sep="_")





datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0.03) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
sOut <- "plot_s0.006-0.02.h0.AFtrace.pdf";

datSimsToPlot <- datSimsToPlot[order(-datSimsToPlot$TopPercentileForBreeding,-datSimsToPlot$GenotypedMarkerPerc, datSimsToPlot$s, -datSimsToPlot$h, datSimsToPlot$s_sd, -datSimsToPlot$h_sd, datSimsToPlot$PopSize), ];

arrColors <- as.character(wes_palette("Zissou1", 6, type = "continuous"))

pdf(file=sOut, width=10, height=5 );


  
  # Save current graphical parameters
  opar <- par(no.readonly = TRUE)
  # Margins of the plot (the first is the bottom margin)
  par(mar = c(4.1, 4.1, 4.1, 15))
  
 

  datSumSim <- NULL;

  
arrIncludeRows <- c();
datPlotDataSum <- NULL;
for(nRow in 1:nrow(datSimsToPlot)) {
  arrFiles <- Sys.glob(datSimsToPlot$OutStem[nRow])
  arrAFTraceFiles <- Sys.glob(datSimsToPlot$OutStem2[nRow])
  
  if (length(arrFiles) != length(arrAFTraceFiles)) {
    cat("Error, simulation ",datSimsToPlot$OutStem[nRow]," do not have matching number of perlocusfitness and AFtrace files\n");
    next;
  }
  
  if (length(arrFiles)<2) {
    cat("Error, simulation ",datSimsToPlot$OutStem[nRow]," does not have enough repetitions.\n");
    next;
  } else {
    cat("INFO: simulation ",datSimsToPlot$OutStem[nRow]," had ", length(arrFiles) ," repetitions.\n");
  }
  arrIncludeRows <- c(arrIncludeRows, nRow);
  datSimReps <- NULL;
  arrS <- c();

  datFitnessAFTrace <- NULL;
  for(nRep in 1:length(arrFiles) ) {
    sFitnessDef <- arrFiles[nRep]
    sAFTrace <- arrAFTraceFiles[nRep]
    datFitnessDef <- read.table(sFitnessDef, header=T) 
    datAFTrace <- read.table(sAFTrace, header=F)
    if (nrow(datAFTrace)<200) {
      cat("Warning, simulation ",datSimsToPlot$OutStem[nRow], " rep ", nRep, " went extinct at gen ",nrow(datAFTrace) ,", skip\n");
      next;
    }
    arrGen0Freq <- colSums( datAFTrace)[ 2:ncol(datAFTrace) ];
    if (length(arrGen0Freq) != nrow(datFitnessDef)) {
      cat("Error, simulation ",datSimsToPlot$OutStem[nRow]," perlocusfitness and AFtrace files do not have same number of loci\n");
      next;
    }
    arrSegSites <- which(arrGen0Freq>0);
    datFitnessDefSeg <- datFitnessDef[arrSegSites, ]
    datAFtraceSeg <- t(datAFTrace[, arrSegSites+1]);
    #hist(datFitnessDefSeg$s, breaks = 20)
    arrS <- c(arrS, datFitnessDefSeg$s);
    datFitnessDefSeg$SimRep <- nRep;
    datFitnessDefSeg <- cbind(datFitnessDefSeg, datAFtraceSeg);
    datFitnessAFTrace <- rbind(datFitnessAFTrace, datFitnessDefSeg);
    # datFitnessDefSeg$s_cat <- cut(datFitnessDefSeg$s, breaks = 6)
    # datFitnessDefSeg <- cbind(datFitnessDefSeg, datAFtraceSeg);
    # datAFTraceMean <- rbind( datAFTraceMean, aggregate(datFitnessDefSeg[, 6:length(datFitnessDefSeg)], by=list(s_cat=datFitnessDefSeg$s_cat), FUN=function(x) {mean(x)} ));
    # datAFTraceSD <- rbind(datAFTraceSD, aggregate(datFitnessDefSeg[, 6:length(datFitnessDefSeg)], by=list(s_cat=datFitnessDefSeg$s_cat), FUN=function(x) {sd(x)} ));
    
  }
  
  datFitnessRelAFTrace <- datFitnessAFTrace[, 1:5]
  datFitnessRelAFTrace$s_cat <- cut(datFitnessRelAFTrace$s, breaks = 6)
  datFitnessRelAFTrace <- cbind(datFitnessRelAFTrace, datFitnessAFTrace[, 6:length(datFitnessAFTrace)] - datFitnessAFTrace[, 6])
  datRelAFTraceMean <- aggregate(datFitnessRelAFTrace[, 7:length(datFitnessRelAFTrace)], by=list(s_cat=datFitnessRelAFTrace$s_cat, simrep=datFitnessRelAFTrace$SimRep), FUN=function(x) {mean(x)} )
  datRelAFTraceMeanOfMean <- aggregate(datRelAFTraceMean[, 3:length(datRelAFTraceMean)], by=list(s_cat=datRelAFTraceMean$s_cat), FUN=function(x) {mean(x)} )
  datRelAFTraceSDOfMean <- aggregate(datRelAFTraceMean[, 3:length(datRelAFTraceMean)], by=list(s_cat=datRelAFTraceMean$s_cat), FUN=function(x) {sd(x)} )
  
  
  datFitnessAbsAFTrace <- datFitnessAFTrace[, 1:5]
  datFitnessAbsAFTrace$s_cat <- cut(datFitnessAbsAFTrace$s, breaks = 6)
  datFitnessAbsAFTrace <- cbind(datFitnessAbsAFTrace, datFitnessAFTrace[, 6:length(datFitnessAFTrace)])
  datAbsAFTraceMean <- aggregate(datFitnessAbsAFTrace[, 7:length(datFitnessAbsAFTrace)], by=list(s_cat=datFitnessAbsAFTrace$s_cat, simrep=datFitnessAbsAFTrace$SimRep), FUN=function(x) {mean(x)} )
  datAbsAFTraceMeanOfMean <- aggregate(datAbsAFTraceMean[, 3:length(datAbsAFTraceMean)], by=list(s_cat=datAbsAFTraceMean$s_cat), FUN=function(x) {mean(x)} )
  datAbsAFTraceSDOfMean <- aggregate(datAbsAFTraceMean[, 3:length(datAbsAFTraceMean)], by=list(s_cat=datAbsAFTraceMean$s_cat), FUN=function(x) {sd(x)} )
  
  hist(arrS, breaks = 10, main="", xlab=paste0("Distribution of simulated s (mean=",datSimsToPlot[nRow, 's'], ", sd=", datSimsToPlot[nRow, 's_sd']) )
  bNewPlot <- T;
  nCat <- 1;
  arrYLims <- c(min(datRelAFTraceMean[,3:length(datRelAFTraceMean)]) , max(datRelAFTraceMean[,3:length(datRelAFTraceMean)]));
  for(sCat in datRelAFTraceMeanOfMean$s_cat) {
    datSimForPlot <- datRelAFTraceMeanOfMean[datRelAFTraceMeanOfMean$s_cat==sCat, ]
    datSimForPlot2 <- datRelAFTraceSDOfMean[datRelAFTraceSDOfMean$s_cat==sCat, ]
    datSimForPlot <- data.frame(gen=0:(ncol(datSimForPlot)-2), delperhap_mean_mean=as.numeric(datSimForPlot[,2:ncol(datSimForPlot)]) , delperhap_mean_sd=as.numeric(datSimForPlot2[,2:ncol(datSimForPlot2)]) )
    
    if (bNewPlot ) {
      plot(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean, type='l', col=paste0(arrColors[nCat], 'ff'), lwd=2, ylim = arrYLims , main="", xlab="Generations", ylab=sPlotVarName)
      bNewPlot <- F;
      datPlotDataSum <- datSimForPlot;
    } else {
      lines(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean,col=paste0(arrColors[nCat], 'ff'), lwd=2);
      datPlotDataSum <- cbind(datPlotDataSum, datSimForPlot);
    }
    
    polygon(c(rev(datSimForPlot$gen), datSimForPlot$gen), c(rev(datSimForPlot$delperhap_mean_mean + datSimForPlot$delperhap_mean_sd*1.96), datSimForPlot$delperhap_mean_mean - datSimForPlot$delperhap_mean_sd*1.96), col = paste0(arrColors[nCat],'11'), border = NA)
    nCat <- nCat + 1;
  }

  bNewPlot <- T;
  nCat <- 1;
  arrYLims <- c(min(datAbsAFTraceMean[,3:length(datAbsAFTraceMean)]) , max(datAbsAFTraceMean[,3:length(datAbsAFTraceMean)]));
  
  for(sCat in datAbsAFTraceMeanOfMean$s_cat) {
    datSimForPlot <- datAbsAFTraceMeanOfMean[datAbsAFTraceMeanOfMean$s_cat==sCat, ]
    datSimForPlot2 <- datAbsAFTraceSDOfMean[datAbsAFTraceSDOfMean$s_cat==sCat, ]
    datSimForPlot <- data.frame(gen=0:(ncol(datSimForPlot)-2), delperhap_mean_mean=as.numeric(datSimForPlot[,2:ncol(datSimForPlot)]) , delperhap_mean_sd=as.numeric(datSimForPlot2[,2:ncol(datSimForPlot2)]) )
    
    if (bNewPlot ) {
      plot(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean, type='l', col=paste0(arrColors[nCat], 'ff'), lwd=2, ylim = arrYLims , main="", xlab="Generations", ylab=sPlotVarName)
      bNewPlot <- F;
      datPlotDataSum <- datSimForPlot;
    } else {
      lines(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean,col=paste0(arrColors[nCat], 'ff'), lwd=2);
      datPlotDataSum <- cbind(datPlotDataSum, datSimForPlot);
    }
    
    polygon(c(rev(datSimForPlot$gen), datSimForPlot$gen), c(rev(datSimForPlot$delperhap_mean_mean + datSimForPlot$delperhap_mean_sd*1.96), datSimForPlot$delperhap_mean_mean - datSimForPlot$delperhap_mean_sd*1.96), col = paste0(arrColors[nCat],'11'), border = NA)
    nCat <- nCat + 1;
  }
  
  #Plot legend

  legend(x = "topright",
         inset = c(-0.45, 0),
         legend = levels(datAbsAFTraceMeanOfMean$s_cat), 
         col = arrColors,
         lwd = 4, xpd = TRUE
  )
  

  
}


par(opar)


dev.off();