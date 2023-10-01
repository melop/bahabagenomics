setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf");
library(ggplot2)
par(family="Arial")
sRefSp <- "BahahaTaipingensis"

bOnlySegregating <- F;
dat <- read.table("consurf_af_sv.rec.bed", header=F, sep="\t")
datPresence <- read.table("../hifiasm_haplotype_gene_variants/counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)
datOrthID2RNAID <- read.table("../hifiasm_bahaha_haplotypes/synorthos.txt", header=T, sep="\t")

if (bOnlySegregating) {
  dat <- dat[rowSums(dat[,8:11])<4,]
  dat <- dat[complete.cases(dat),];
} else {
  dat <- dat[complete.cases(dat),];
}

datPresence$UniqID <- paste(datPresence$pgChr, datPresence$pgOrd, datPresence$pgID, sep="_")
datOrthID2RNAID$UniqID <- paste(datOrthID2RNAID$pgChr, datOrthID2RNAID$pgOrd, datOrthID2RNAID$pgID, sep="_")
datPresence <- merge(datPresence, datOrthID2RNAID[, c("UniqID",  sRefSp)], by="UniqID", all.x = T, all.y=F);
colnames(datPresence)[ncol(datPresence)] <- 'mRNAID'

# dat <- dat[ (dat$V8 + dat$V9 + dat$V10 + dat$V11 != 4), ] # exclude all sites that are completely fixed 
# dat <- dat[!is.na(dat$V1), "V8"];

arrHighImpactVarTypes <- c("frameshift_variant", "start_lost" , "stop_gained", "stop_lost", "exon_loss_variant" ,            
                           "splice_donor_variant", "splice_acceptor_variant"  ,      
                             "gene_fusion"  , "ablation"   );


datByGene <- aggregate(dat$V15, by=list(gene=dat$V4), FUN=mean)
arrHighGenes <-  unique(dat[dat$V6 == 'HIGH' , 'V4']);
arrLowGenes <-  unique(dat[dat$V6 == 'LOW' , 'V4']);

arrHighRecRate <- datByGene[datByGene$gene %in% arrHighGenes, 'x']
arrLowRecRate <- datByGene[(!(datByGene$gene %in% arrHighGenes)) & (datByGene$gene %in% arrLowGenes), 'x']

median(arrHighRecRate);
median(arrLowRecRate);
wilcox.test(arrHighRecRate, arrLowRecRate)

mean(log10(arrHighRecRate) ) ;
mean(log10(arrLowRecRate) );
t.test(log10(arrHighRecRate), log10(arrLowRecRate))

#check SNP mutations that changes single AA:

datConsurfRec <- read.table("consurf_af.rec.bed", header=F, sep="\t")
datConsurfRecPergene <- aggregate(datConsurfRec[,c('V5', 'V13')], by=list(gene=datConsurfRec$V4), FUN=mean)
datConsurfRecPergene <- datConsurfRecPergene[!(datConsurfRecPergene$gene %in% arrHighGenes), ]; #remove genes with high impact SVs
summary(lm(log10(datConsurfRecPergene$V13) ~ datConsurfRecPergene$V5));

hist(datConsurfRec[datConsurfRec$V6==1 ,'V5' ]);
hist(datConsurfRec[datConsurfRec$V7==1 ,'V5' ]);
hist(datConsurfRec[datConsurfRec$V8==1 ,'V5' ]);
hist(datConsurfRec[datConsurfRec$V9==1 ,'V5' ]);

#do it haplotype by haplotype:

datPlot <- NULL;
datVarTypePlot <- NULL;
for(nHap in 1:4) {
  cat("========== HAP ", nHap , "=============\n");
  datByGene <- dat[dat[,7+nHap]==1,]
  datByGene <- aggregate(datByGene$V15, by=list(gene=datByGene$V4), FUN=mean)
  arrHighGenes <-  unique(dat[dat$V6 == 'HIGH' , 'V4']);
  arrLowGenes <-  unique(dat[dat$V6 == 'LOW' , 'V4']);
  
  datHigh <- dat[dat$V6 == 'HIGH' & dat[,7+nHap]==1, c('V4' ,'V5', 'V15') ]
  datHigh$hap <- nHap;
  for(sType in arrHighImpactVarTypes) {
    nTypeCount <-  length(grep(sType, datHigh$V5, fixed=T));
    if (sType == "ablation") {#add gene deletion counts
      arrMoreMissingGenes <- datPresence[datPresence[,5+nHap]==0, 'mRNAID'];
      arrHighGenes <- c(arrHighGenes, arrMoreMissingGenes);
      nMoreCount <-  sum(as.integer( !(arrMoreMissingGenes %in% datHigh$V4) ));
      nTypeCount <- nTypeCount + nMoreCount
    }
    datVarTypePlot <- rbind(datVarTypePlot, data.frame(hap = nHap, vartype=sType, count=nTypeCount) );
  }
  
  
  arrHighRecRate <- datByGene[datByGene$gene %in% arrHighGenes, 'x']
  arrLowRecRate <- datByGene[(!(datByGene$gene %in% arrHighGenes)) & (datByGene$gene %in% arrLowGenes), 'x']
  
  cat('Genes with high impact median rec rate ', median(arrHighRecRate),"\n")
  cat('Genes with only low impact median rec rate ',median(arrLowRecRate),"\n");
  print(wilcox.test(arrHighRecRate, arrLowRecRate))
  
  cat('Genes with high impact mean log10 rec rate ', mean(log10(arrHighRecRate)), "\n");
  cat('Genes with only low impact mean log10 rec rate ', mean(log10(arrLowRecRate) ) ,"\n");
  print(t.test(log10(arrHighRecRate), log10(arrLowRecRate) ))
  
  datPlot <- rbind(datPlot, data.frame(hap=nHap, impact="high", generecrate=arrHighRecRate, stringsAsFactors = F) );
  datPlot <- rbind(datPlot, data.frame(hap=nHap, impact="low", generecrate=arrLowRecRate, stringsAsFactors = F) );
  
  datConsurfRec <- read.table("consurf_af.rec.bed", header=F, sep="\t")
  # datConsurfRec$fixedder <- (rowSums(datConsurfRec[, (5+1:4)])>1) #high freq
  # datConsurfRec <- datConsurfRec[!datConsurfRec$fixedder, ]; #remove high freq
  # datConsurfRec <- datConsurfRec[complete.cases(datConsurfRec),];
  datConsurfRec <- datConsurfRec[datConsurfRec[,5+nHap]==1,]
  datConsurfRecPergene <- aggregate(datConsurfRec[,c('V5', 'V13')], by=list(gene=datConsurfRec$V4), FUN=mean)
  datConsurfRecPergene <- datConsurfRecPergene[!(datConsurfRecPergene$gene %in% arrHighGenes), ]; #remove genes with high impact SVs
  
  cat("Per gene log10(rho) ~ mean_consurf\n");
  print(summary(lm(log10(datConsurfRecPergene$V13) ~ datConsurfRecPergene$V5)));
  cat("Per site log10(rho) ~ consurf\n");
  print(summary(lm(log10(datConsurfRec$V13) ~ datConsurfRec$V5)));
  
  datConsurfRec$rec <- datConsurfRec$V13/100;
  datConsurfHighCountRecPergene <- aggregate(datConsurfRec[,c('V5', 'rec')], by=list(gene=datConsurfRec$V4), FUN=function(x) { if(x[1]>=1) {sum(as.integer(x>=9))/length(x)} else {mean(x)} })
  datConsurfHighCountRecPergene <- datConsurfHighCountRecPergene[!(datConsurfHighCountRecPergene$gene %in% arrHighGenes), ]; #remove genes with high impact SVs
  
  cat("Per gene log10(rho) ~ high consurf site count percentage\n");
  print(summary(lm(log10(datConsurfHighCountRecPergene$rec) ~ datConsurfHighCountRecPergene$V5)));
  
  
}

datPlot$hap <- as.factor(datPlot$hap)
datPlot$impact <- as.factor(datPlot$impact)

datPlot$hapL <- factor(datPlot$hap, levels = 1:4, labels = c('F0', 'F1', 'M0', 'M1') )
#boxplot(log10(generecrate) ~ impact * hap , data=datPlot, pch=16, col=c('#3B9AB2' , '#F21A00') );

ggplot(datPlot, aes(x=hapL, y=log10(generecrate), fill=impact)) + 
  geom_boxplot()+
  theme_minimal() + 
  scale_fill_manual('High Impact Variant Type', values=rev(as.character(wes_palette("Zissou1", 2, type = "continuous"))) )

if (bOnlySegregating) {
  ggsave("impacttype_vs_recrate.segregating.pdf", width=6, height=6)
} else {
  ggsave("impacttype_vs_recrate.pdf", width=6, height=6)
}
summary(lm(log10(generecrate) ~ impact + hap, data=datPlot));

ggplot(datVarTypePlot, aes(fill=vartype, y=count, x=hap)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Count', title='High impact variants per haplotype') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(arrHighImpactVarTypes), type = "continuous")) )

if (bOnlySegregating) {
  ggsave("impacttypes_per_hap.segregating.pdf", width=6, height=6)
} else {
  ggsave("impacttypes_per_hap.pdf", width=6, height=6)
}
summary(lm(datVarTypePlot$count ~ datVarTypePlot$hap + datVarTypePlot$vartype))

