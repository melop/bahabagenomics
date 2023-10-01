setwd("/data/projects/rcui/bahaha_assembly/hyphy/relax_genespace_nooutgroup/rho");
# library(ggplot2)
#datRecRate <- read.table("/data/projects/rcui/bahaha_assembly/ismc/20ind/20ind/gene_recrate.aggregated.txt", header = T, stringsAsFactors = F)

datRecRate <- read.table("/data/projects/rcui/bahaha_assembly/ismc/20ind/20ind/gene_recrate.aggregated.excludeROH.txt", header = T, stringsAsFactors = F)

datRelax <- read.table("../sum_BTP.txt", header=F, sep="\t", fill = T, quote = "")
datRelax2 <- read.table("../../relax_genespace/sum_BTP.txt", header=F, sep="\t", fill = T, quote = "")

datRelaxM <- merge(datRelax[, c(1,5,9)], datRelax2[, c(1,5,9)], by="V1", all.x=T, all.y=T)

datRelaxM$V5 <- 1;
datRelaxM$V9 <- 0;

for(i in 1:nrow(datRelaxM)) {
  if (is.na(datRelaxM[i , 'V5.x'])) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.y'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.y'];
  } else if (is.na(datRelaxM[i , 'V5.y'])) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.x'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.x'];
  } else if (datRelaxM[i , 'V5.y'] < datRelaxM[i , 'V5.x'] ) {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.y'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.y'];   
  } else {
    datRelaxM[ i , 'V5'] <- datRelaxM[i , 'V5.x'];
    datRelaxM[ i , 'V9'] <- datRelaxM[i , 'V9.x'];    
  }
}

datRelax <- datRelaxM[, c(1, 6, 7)];

datRelax <- merge(datRelax, datRecRate, by.x='V1', by.y="OrthoID")

datRelax$V9[datRelax$V9<1/50] <- 1/50; #the max k is 50, so the minimal k should be set to 1/50 to be fair
datRelax$logK <- log10(datRelax$V9)
summary(lm(datRelax$logK ~ log10(datRelax$recrate)) )

arrRelaxedRec <- datRelax[datRelax$V5<0.001 & datRelax$V9 <1, "recrate"]
arrIntenseRec <- datRelax[datRelax$V5<0.001 & datRelax$V9 >1, "recrate"]
arrNonRelaxedRec <- datRelax[datRelax$V5>=0.001 , "recrate"]
wilcox.test(arrRelaxedRec, arrNonRelaxedRec, paired = F)
wilcox.test(arrRelaxedRec, arrIntenseRec, paired = F)
wilcox.test(arrIntenseRec, arrNonRelaxedRec, paired = F)
median(arrRelaxedRec)
median(arrIntenseRec)
median(arrNonRelaxedRec)

mean(arrRelaxedRec)
mean(arrIntenseRec)
mean(arrNonRelaxedRec)

datPlot <- data.frame(recrate=c(arrRelaxedRec,arrIntenseRec,arrNonRelaxedRec ), relaxtype=c(rep("relax", length(arrRelaxedRec)), rep("intense", length(arrIntenseRec) ),rep("other", length(arrNonRelaxedRec) ) ) )
datPlot$relaxtype <- as.factor(datPlot$relaxtype)

pdf(file="relax_rec.pdf", width=3, height = 6)
boxplot(datPlot$recrate ~ datPlot$relaxtype, ylim=c(0,0.4), col=c( "#F21A00", "#D1C74C", "#3B9AB2" ), pch=16 )


data_means <- aggregate(datPlot$recrate,                  
                        list(datPlot$relaxtype),
                        mean)

points(x = 1:nrow(data_means),                           
       y = data_means$x,
       col = c("blue", "red", "red"),
       pch = 16)

# p<-ggplot(datPlot, aes(x = relaxtype, y = recrate)) 
# p+geom_violin()

dev.off();
