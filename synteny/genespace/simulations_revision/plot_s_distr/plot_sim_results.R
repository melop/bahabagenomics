setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/simulations_revision/plot_s_distr");
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


arrPlotVars <- c("delperhap_mean", "selection_killedind", "fitness_beforesel_mean", "fitness_beforesel_sd");
arrPlotVarNames <- c("Deleterious mutations per haplotype", 
                     "Proportion of individuals killed due to viability selection", 
                     "Mean absolute fitness", "Stddev of absolute fitness");
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
arrDefs <- Sys.glob("simdef*.txt")
#arrDefs <- c("simdef.txt", "simdef2.txt")

for(sDef in arrDefs) {
  datSimDefs <- rbind(datSimDefs, read.table(sDef, header=T, sep="\t") )
}

datSimDefs$OutStem <- paste("../*/out/simret","gen",datSimDefs$TotalGen,
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
                  'seed', '[0-9]*.txt', sep="_")


# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002,0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.5, 0.25) & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_main1.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1,0.5, 0.25, 0.1, 0.05), ]
# sOut <- "plot_main2.pdf";

#datSimsToPlot <- datSimDefs[datSimDefs$TopPercentileForBreeding %in% c(1,0.25),]

#datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100,1000) & datSimDefs$h == 0 & datSimDefs$s == 0.006 & datSimDefs$s_sd == 0.006/2 & datSimDefs$TopPercentileForBreeding %in% c(1) & datSimDefs$GenotypedMarkerPerc %in% c(1),]

#datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100,1000) & datSimDefs$h == 0 & datSimDefs$h_sd == 0 & datSimDefs$s_sd == 0 & datSimDefs$s %in% c(0.001, 0.002, 0.003, 0.004, 0.005) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001,0.002,0.003,0.004,0.005,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0_hsd_0.25_no_artf.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001,0.002,0.003,0.004,0.005,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_no_artf.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_ssd_0_var_no_artf.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001,0.002,0.003,0.004,0.005,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0_hsd_0.25_no_artf_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001,0.002,0.003,0.004,0.005,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_no_artf_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding == 1 & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_ssd_0_var_no_artf_N1000.pdf";


# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.5, 0.25) & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0_hsd_0.25_no_artf_K0.5.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.5, 0.25) & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_no_artf_K0.5.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.5, 0.25) & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0_hsd_0.25_no_artf_K0.5_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002,0.004,0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.5, 0.25) & datSimDefs$GenotypedMarkerPerc==1, ]
# sOut <- "plot_purifying_h0.5_hsd_0.25_no_artf_K0.5_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.001_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.002_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.003_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.004) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.004_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.005) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.005_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.006_h0_hsd_0.25_no_artf_K0.25_varyingM.pdf";


# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.001_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.002_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.004) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.004_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.005) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.005_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.006_h0.5_hsd_0.25_no_artf_K0.25_varyingM.pdf";



# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.001_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.002) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05), ]
# sOut <- "plot_purifying_s0.002_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.004) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.004_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.005) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.005_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(1000) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.5, 0.25, 0.1 ,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_purifying_s0.006_h0.5_hsd_0.25_no_artf_K0.25_varyingM_N1000.pdf";


# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(1) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K1_varyingSDDel_N100.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(0.5) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.5_varyingSDDel_N100.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c( 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.25_varyingSDDel_N100.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(100) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(1) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K1_varyingSDDel_N100_N0100.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(100) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(0.5) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.5_varyingSDDel_N100_N0100.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(100) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.003) & datSimDefs$TopPercentileForBreeding %in% c(0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1) & datSimDefs$SDDel>0, ]
# sOut <- "plot_purifying_s0.003_h0.5_hsd_0.25_no_artf_K0.25_varyingSDDel_N100_N0100.pdf";

#datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.006) & datSimDefs$TopPercentileForBreeding %in% c(1,0.5,0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1) , ]
#datSimsToPlot <- rbind(datSimsToPlot, datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h >=0 & datSimDefs$h_sd==0 & datSimDefs$s_sd == 0 & datSimDefs$s %in% c(0) & datSimDefs$TopPercentileForBreeding %in% c(1,0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1,0.05) , ])

#datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>0 & datSimDefs$s_sd > 0 & datSimDefs$s %in% c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006) & datSimDefs$TopPercentileForBreeding %in% c(1) & datSimDefs$GenotypedMarkerPerc %in% c(1) , ]

#datSimsToPlot <- datSimDefs[datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.005) & datSimDefs$TopPercentileForBreeding %in% c(1 , 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1,0.5, 0.05), ]
#sOut <- "plot1.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.006) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_neteffect_s0.006.pdf";
# 
# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.005) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_neteffect_s0.005.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.004) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_neteffect_s0.004.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.003) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_neteffect_s0.003.pdf";

# datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0.5 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0, 0.002) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
# sOut <- "plot_neteffect_s0.002.pdf";

datSimsToPlot <- datSimDefs[datSimDefs$SeedingPopSize %in% c(10) & datSimDefs$PopSize %in% c(100) & datSimDefs$h == 0 & datSimDefs$h_sd>=0 & datSimDefs$s_sd >= 0 & datSimDefs$s %in% c(0.006, 0.01, 0.02) & datSimDefs$TopPercentileForBreeding %in% c(1, 0.25) & datSimDefs$GenotypedMarkerPerc %in% c(1, 0.25,0.05) & datSimDefs$SDDel==-1, ]
sOut <- "plot_s0.006-0.02.h0.pdf";

datSimsToPlot <- datSimsToPlot[order(-datSimsToPlot$TopPercentileForBreeding,-datSimsToPlot$GenotypedMarkerPerc, datSimsToPlot$s, -datSimsToPlot$h, datSimsToPlot$s_sd, -datSimsToPlot$h_sd, datSimsToPlot$PopSize), ];

arrColors <- as.character(wes_palette("Zissou1", nrow(datSimsToPlot), type = "continuous"))

pdf(file=sOut, width=10, height=5 );

for(nVar in 1:length(arrPlotVars)) {
  
  # Save current graphical parameters
  opar <- par(no.readonly = TRUE)
  # Margins of the plot (the first is the bottom margin)
  par(mar = c(4.1, 4.1, 4.1, 15))
  
  
  ylim <- c(0,arrYLimMax[nVar]);
  sPlotVar <- arrPlotVars[nVar];
  sPlotVarName <- arrPlotVarNames[nVar];
  datSumSim <- NULL;
  bNewPlot <- T;
  
arrIncludeRows <- c();
datPlotDataSum <- NULL;
for(nRow in 1:nrow(datSimsToPlot)) {
  arrFiles <- Sys.glob(datSimsToPlot$OutStem[nRow])
  arrFiles <- arrFiles[!grepl('perlocusfitness', arrFiles)]
  if (length(arrFiles)<2) {
    cat("Error, simulation ",datSimsToPlot$OutStem[nRow]," does not have enough repetitions.\n");
    next;
  } else {
    cat("INFO: simulation ",datSimsToPlot$OutStem[nRow]," had ", length(arrFiles) ," repetitions.\n");
  }
  arrIncludeRows <- c(arrIncludeRows, nRow);
  datSimReps <- NULL;
  for(sF in arrFiles) {
    datSim <- read.table(sF, header=T) 
    if (is.null(datSimReps)) {
      datSimReps <- datSim[, c('gen', sPlotVar)];
    } else {
      datSimReps <- merge(datSimReps, datSim[, c('gen', sPlotVar)], by="gen", all.x=T, all.y=F)
    }
  }
  
  datSimForPlot <- data.frame(gen=datSimReps$gen, delperhap_mean_mean=rowMeans(datSimReps[,2:ncol(datSimReps)]) , delperhap_mean_sd=apply(datSimReps[,2:ncol(datSimReps)], 1, sd) )
  datSimForPlot <- datSimForPlot[complete.cases(datSimForPlot) & datSimForPlot$gen>0, ];
  if (bNewPlot ) {
    plot(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean, type='l', col=paste0(arrColors[nRow], 'ff'), lwd=2, xlim=xlim, ylim = ylim , main="", xlab="Generations", ylab=sPlotVarName)
    bNewPlot <- F;
    datPlotDataSum <- datSimForPlot;
  } else {
    lines(datSimForPlot$gen, datSimForPlot$delperhap_mean_mean,col=paste0(arrColors[nRow], 'ff'), lwd=2);
    datPlotDataSum <- cbind(datPlotDataSum, datSimForPlot);
  }
  
  polygon(c(rev(datSimForPlot$gen), datSimForPlot$gen), c(rev(datSimForPlot$delperhap_mean_mean + datSimForPlot$delperhap_mean_sd*1.96), datSimForPlot$delperhap_mean_mean - datSimForPlot$delperhap_mean_sd*1.96), col = paste0(arrColors[nRow],'55'), border = NA)
  
  
}

#Plot legend
#Figure out what parameters are shared between sims
datSimsToPlot <- datSimsToPlot[arrIncludeRows,];
arrColors <- arrColors[arrIncludeRows]
arrParamIdentical <- apply(datSimsToPlot[,3:15], 2, function(x) {length(unique(x))==1})

arrSharedParamNames <- names(arrParamIdentical)[arrParamIdentical]
arrSharedParamValues <- as.character(datSimsToPlot[1,arrSharedParamNames])
arrRemoveDefaults <- arrSharedParamValues!="-1";
arrSharedParamValues <- arrSharedParamValues[arrRemoveDefaults] 
arrSharedParamNames <- arrSharedParamNames[arrRemoveDefaults]
arrSharedLabels <- paste(arrParamNameMap[arrSharedParamNames], arrSharedParamValues, sep=" = ")

arrDiffParamNames <- names(arrParamIdentical)[!arrParamIdentical]
datSimDiffParams <- datSimsToPlot[,arrDiffParamNames, drop=F]
arrDiffParamLabels <- apply(datSimDiffParams, 1, function(x) {paste(arrParamNameMap[arrDiffParamNames], x, sep="=", collapse = ', ')})
# Add first legend
legend(x = "topright",
       inset = c(-0.45, 0),
       legend = arrDiffParamLabels, 
       col = arrColors,
       lwd = 4, xpd = TRUE
)

legend(x = "bottomright",
       inset = c(-0.45, 0),
       legend = arrSharedLabels, 
       col = 'black',
       seg.len=0, xpd = TRUE
)
par(opar)
}

dev.off();