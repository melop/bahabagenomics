setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants_consurf");
library(ggplot2)
library(wesanderson)

arrHaps <- c('F0', 'F1', 'M0', 'M1');
bExcludeFixed <- T;


arrHighImpactVarTypes <- c("ablation", "frameshift_variant", "5_prime_UTR_premature_start_codon_gain_variant", "start_lost" , "stop_gained", "stop_lost", "exon_loss_variant" ,            
                           "splice_donor_variant", "splice_acceptor_variant"  );
names(arrHighImpactVarTypes) <- gsub('_variant','', arrHighImpactVarTypes);
names(arrHighImpactVarTypes) <- gsub('_',' ', names(arrHighImpactVarTypes) );
names(arrHighImpactVarTypes)[3] <- "start gain"


datVarTypePlot2 <- read.table("consurf.SV.polarized.out.txt", header=F, sep="\t")

if (bExcludeFixed) {
  datVarTypePlot2 <- datVarTypePlot2[rowSums(datVarTypePlot2[,(1:length(arrHaps)+8) ] ) < length(arrHaps), ];
}

datVarTypePlot2 <- datVarTypePlot2[datVarTypePlot2$V7 == 'HIGH', ]
datVarTypePlotHighImpact <- NULL;

for(nHap in 1:length(arrHaps)) {
  
  sHap <- arrHaps[nHap];
  datVars <- datVarTypePlot2[, c(6, 8+nHap)];
  for(sType in arrHighImpactVarTypes) {
    nCount <- sum(grepl(sType, datVars[,1]) & (datVars[,2]==1) , na.rm = T);
    if (nCount > 0) {
      datRow <- data.frame(hap=sHap, vartype=sType, count=nCount);
      datVarTypePlotHighImpact <- rbind(datVarTypePlotHighImpact, datRow);
    }
  }
  
}

for(i in 1:length(arrHighImpactVarTypes)) {
  datVarTypePlotHighImpact$vartype <- str_replace_all(as.character(datVarTypePlotHighImpact$vartype), as.character(arrHighImpactVarTypes)[i] , names(arrHighImpactVarTypes)[i] )
}


ggplot(datVarTypePlotHighImpact, aes(fill=vartype, y=count, x=hap)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Count', title='High impact variants per haplotype') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(unique(datVarTypePlotHighImpact$vartype)), type = "continuous")) )
if (bExcludeFixed) {
  ggsave("impacttypes_highimpact_per_hap.directaln.exclude.fixedsites.pdf", width=6, height=6)
  
} else {
  ggsave("impacttypes_highimpact_per_hap.directaln.pdf", width=6, height=6)
}


mat <- matrix(c(datVarTypePlotHighImpact$count[1:16]), ncol=4, byrow=TRUE )
chisq.test(t(mat) )


mat <- matrix(c(datVarTypePlotHighImpact$count[1:4], datVarTypePlotHighImpact$count[17], 
                datVarTypePlotHighImpact$count[5:8], datVarTypePlotHighImpact$count[18],
                datVarTypePlotHighImpact$count[9:12], datVarTypePlotHighImpact$count[19],
                datVarTypePlotHighImpact$count[13:16], datVarTypePlotHighImpact$count[20]
                ), ncol=5, byrow=TRUE )
chisq.test(t(mat) )

