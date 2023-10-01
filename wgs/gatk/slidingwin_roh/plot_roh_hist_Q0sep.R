setwd("/data/projects/rcui/bahaha_assembly/wgs/gatk/slidingwin_roh");
library(chromPlot)
library(RColorBrewer)
library(wesanderson)
library(svglite)

#pdf(file = "fig4.roh.BTP_Q0.revised.pdf",width=10,height=8);
#pdf(file = "fig4.roh.parentalpop.revised.pdf",width=10,height=8);
svglite(file = "hist_roh_q0_sep.svg",width=6,height=5);
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


arrHomFiles <- Sys.glob("*.hom.mid.win.tsv");
#remove M06 since it is same individual as Q1
arrHomFiles <- arrHomFiles[!grepl("BTP_Q1", arrHomFiles)]
#arrHomFiles <- arrHomFiles[grepl("BTP_Q0", arrHomFiles)]
#arrHomFiles <- arrHomFiles[!grepl("BTP_Q0", arrHomFiles)]

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


datROHSegments$len <- datROHSegments$homend - datROHSegments$homstart + 1
hist(log10(datROHSegments$len[datROHSegments$Group!="Q0"]), xlim=c(4,7), breaks=20 , col="#E7C21C88", freq = F, main="ROH lengths", xlab="log10(bp)")
hist(log10(datROHSegments$len[datROHSegments$Group=="Q0"]), xlim=c(4,7), breaks=20 , col="#1CC2E788", freq = F, add=T)

dev.off();
