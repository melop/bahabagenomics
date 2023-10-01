setwd("~/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/rho")
library(vioplot)
library(stringr)
library(ggplot2)
library(wesanderson)
arrBreaks <- c(seq(0.0, 0.8, 0.2) , 0.99, 1);
AF_low <- 0.3
AF_high <- 0.7
dat <- read.table("consurf_af_sv.rec.bed", header=F, sep="\t")
dat <- dat[ , c(1:3, 6:ncol(dat) ) ]
dat <- dat[complete.cases(dat),];

arrHighImpactVarTypes <- c("frameshift_variant", "start_lost" , "stop_gained", "stop_lost", "exon_loss_variant" ,            
                           "splice_donor_variant", "splice_acceptor_variant"  ,      
                           "gene_fusion"  , "ablation" ,"synonymous_variant"  );


dat$FIS <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; nExpHet<- 2*p*(1-p); nObsHet <- length(x[x==1])/length(x) ; if (nExpHet !=0) {FIS <- (nExpHet -nObsHet)/nExpHet; return(FIS)} else {return (0) }  })
dat$HomozygoteBias <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; nMajorHom<-2; nExpHom<- p*p; nObsHom <- length(x[x==nMajorHom])/length(x) ; if (nExpHom !=0) {FIS <- (nExpHom -nObsHom)/nExpHom; return(FIS)} else {return (0) }  })
dat$AF <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; return(p) })
dat <- dat[dat$AF>0, ];
dat$AllGenoIdentical <- apply((dat[, 10:29]), 1, FUN= function(x) {return(length(unique(x))==1 ) })
dat$AFBin <- cut(dat$AF, breaks = arrBreaks);

#look at frameshifting :

datFS <- dat[grepl("frameshift", dat$V9, fixed = T), ];
datFS$mutation_order <- paste(datFS$V6, datFS$V7, sep = ">")
datFS$is_insertion <- str_length(datFS$V6) < str_length(datFS$V7)
sum(datFS$is_insertion == T)
sum(datFS$is_insertion == F)
sum(str_length(datFS$V6) == 1 & str_length(datFS$V7) == 2)
sum(str_length(datFS$V6) == 2 & str_length(datFS$V7) == 1)
arrMutationTypeCounts <- table(datFS$mutation_order);
arrMutationTypeCounts <- arrMutationTypeCounts[order(-arrMutationTypeCounts)];

datInsertion1bp <- datFS[str_length(datFS$V6) == 1 & str_length(datFS$V7) == 2, ]
datDeletion1bp <- datFS[str_length(datFS$V6) == 2 & str_length(datFS$V7) == 1, ]

datMatInsertion1bp <- data.frame(A = c(0,0,0,0), T = c(0,0,0,0), C = c(0,0,0,0), G = c(0,0,0,0) );
rownames(datMatInsertion1bp) <- c('A', 'T', 'C', 'G')
datMatDeletion1bp <- datMatInsertion1bp

for(sB1 in rownames(datMatInsertion1bp)) {
  for(sB2 in rownames(datMatInsertion1bp)) {
    datMatInsertion1bp[sB2, sB1] <- sum(datInsertion1bp$V6 == sB1 & substring(datInsertion1bp$V7, 2) == sB2 );
    datMatDeletion1bp[sB2, sB1] <- sum(datDeletion1bp$V7 == sB1 & substring(datDeletion1bp$V6, 2) == sB2 );
  }
}

write.table(datMatInsertion1bp, file="mat_insertion1bp.txt" , sep = "\t", col.names = T, row.names = T);
write.table(datMatDeletion1bp, file="mat_deletion1bp.txt" , sep = "\t", col.names = T, row.names = T);

#plot variant types by AF:



data_summary_vartypes <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    arrRet <- c();
    for(sTerm in arrHighImpactVarTypes) {
      arrRet <- c(arrRet,  length(grep(sTerm, x[[col]], fixed=T)) )
    }
    names(arrRet) <- arrHighImpactVarTypes;
    return(arrRet);
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  #data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

dsum <- data_summary_vartypes(dat, 'V9'  , c('AFBin') );
dsum <- dsum[complete.cases(dsum),]

arrSum <- rowSums(dsum[,2:(ncol(dsum)-1)])

datSumPlot <- NULL;
for(sTerm in arrHighImpactVarTypes) {
  if (sTerm == "synonymous_variant") {
    next;
  }
  datAdd <- dsum[, c('AFBin', sTerm)];
  datAdd$vartype <- sTerm;
  datAdd$perc <- datAdd[, sTerm] / arrSum;
  datAdd$count <- datAdd[, sTerm];
  datAdd$div_by_syn <- datAdd[, sTerm] / dsum$synonymous_variant;
  datSumPlot <- rbind(datSumPlot, datAdd[ , c("AFBin", "vartype", 'perc', 'count', 'div_by_syn')] );
}

ggplot(datSumPlot, aes(fill=vartype, y=perc, x=AFBin)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Proportion of variants', title='High impact variants at different AF') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(arrHighImpactVarTypes), type = "continuous")) )

ggsave("impacttypes_percent_AF.pdf", width=7, height=3.5)

# ggplot(datSumPlot, aes(fill=vartype, y=div_by_syn, x=AFBin)) + 
#   geom_bar(position='stack', stat='identity') +
#   theme_minimal() + 
#   labs(x='Team', y='Variant count / syn count', title='High impact variants at different AF') +
#   theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
#   scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(arrHighImpactVarTypes), type = "continuous")) )
# 
# ggsave("impacttypes_div_by_syn_AF.pdf", width=7, height=3.5)


ggplot(datSumPlot, aes(fill=vartype, y=count, x=AFBin)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Count of variants', title='High impact variants at different AF') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(arrHighImpactVarTypes), type = "continuous")) )
ggsave("impacttypes_count_AF.pdf", width=7, height=3.5)

datPlot <- dat[dat$AF > AF_low & dat$AF < AF_high & (!dat$AllGenoIdentical), ];

datPlot$impactlevel <- factor(datPlot$V10, c('LOW', 'MODERATE', 'HIGH') );
pdf(file = "FIS.BAH.by.impactlevel.pdf", width=5,height=4)
vioplot(datPlot$FIS ~ datPlot$impactlevel, col=as.character(wes_palette("Zissou1", 3, type = "continuous")) )
vioplot(datPlot$HomozygoteBias ~ datPlot$impactlevel, col=as.character(wes_palette("Zissou1", 3, type = "continuous")))
dev.off();

summary(lm(datPlot$FIS ~ datPlot$impactlevel))
summary(lm(datPlot$HomozygoteBias ~ datPlot$impactlevel))

dat$impact_level <- factor(dat$V10, c('LOW', 'MODERATE', 'HIGH'));
summary(lm(dat$AF ~ dat$impact_level))

summary(lm(dat$AF ~ dat$V9))

arrHighFIS <- dat[dat$V10=="HIGH" & dat$AF > AF_low & dat$AF < AF_high & (!dat$AllGenoIdentical), 'FIS'];

arrHighGenes <- unique(dat[dat$V10 == "HIGH", 'V8' ])
arrLowFIS <- dat[(!(dat$V8 %in% arrHighGenes)) & dat$V10 == "LOW" & dat$AF > AF_low & dat$AF < AF_high & (!dat$AllGenoIdentical), 'FIS']

median(arrHighFIS);
median(arrLowFIS);
wilcox.test(arrHighFIS, arrLowFIS)


mean(arrHighFIS);
mean(arrLowFIS);
t.test(arrHighFIS, arrLowFIS)

hist(arrHighFIS,breaks=20)
hist(arrLowFIS, breaks=20)

hist(dat[dat$V10 == "HIGH" & dat$AF > AF_low & dat$AF < AF_high, "AF"]);
hist(dat[(!(dat$V8 %in% arrHighGenes)) & dat$V10 == "LOW" & dat$AF > AF_low & dat$AF < AF_high, "AF"]);

#conclusion: High impact variants are more likely to have excessive heterozygosity, suggesting homozygous deleterious genotypes have lower fitness (inviable probably)


arrHighHomBias <- dat[dat$V10=="HIGH" & dat$AF > AF_low & dat$AF < AF_high & (!dat$AllGenoIdentical), 'HomozygoteBias'];

arrHighGenes <- unique(dat[dat$V10 == "HIGH", 'V8' ])
arrLowHomBias <- dat[(!(dat$V8 %in% arrHighGenes)) & dat$V10 == "LOW" & dat$AF > AF_low & dat$AF < AF_high & (!dat$AllGenoIdentical), 'HomozygoteBias']

median(arrHighHomBias);
median(arrLowHomBias);
wilcox.test(arrHighHomBias, arrLowHomBias)


mean(arrHighHomBias);
mean(arrLowHomBias);
t.test(arrHighHomBias, arrLowHomBias)

#conclusion: high impact variants are more likely to be biased towards homozygous ancestral genotypes, than homozygous derived, than low impact variants.


datByGene <- aggregate(dat$V35, by=list(gene=dat$V8), FUN=mean)
arrHighGenes <-  unique(dat[dat$V10 == 'HIGH' , 'V8']);
arrLowGenes <-  unique(dat[dat$V10 == 'LOW' , 'V8']);

arrHighRecRate <- datByGene[datByGene$gene %in% arrHighGenes, 'x']
arrLowRecRate <- datByGene[(!(datByGene$gene %in% arrHighGenes)) & (datByGene$gene %in% arrLowGenes), 'x']

median(arrHighRecRate);
median(arrLowRecRate);
wilcox.test(arrHighRecRate, arrLowRecRate)

mean(log10(arrHighRecRate) ) ;
mean(log10(arrLowRecRate) );
t.test(log10(arrHighRecRate), log10(arrLowRecRate))


datRecPlot <- dat
ggplot(datRecPlot, aes(x=AFBin, y=log10(V35), fill=impact_level)) + 
  geom_boxplot()+
  theme_minimal() + 
  scale_fill_manual('High Impact Variant Type', values=(as.character(wes_palette("Zissou1", 3, type = "continuous"))) )
ggsave("AF_impacttype_vs_recrate.pdf", width=8, height=6)
summary(lm(log10(V35) ~ impact_level + AFBin , data = datRecPlot))
#rec rate is lower in genes with high impact variants

