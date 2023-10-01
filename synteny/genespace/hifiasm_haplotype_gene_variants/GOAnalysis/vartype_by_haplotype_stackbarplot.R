setwd("/public3/group_crf/home/cuirf/bahaha_assembly/synteny/genespace/hifiasm_haplotype_gene_variants/GOanalysis");
library(ggplot2)
library(wesanderson)

datPresence <- read.table("../counts_deleterious_ret.presence.txt", header = T, stringsAsFactors = F)

datPCounts <- datPresence[,5:8]
datPCounts[datPCounts==0] <- 10
datPCounts[datPCounts!=10] <- 0
datPCounts[datPCounts==10] <- 1

arrGeneLossCounts <- colSums(datPCounts)

arrHighImpactVarTypes <- c("ablation", "frameshift_variant",  "start_lost" , "stop_gained", "stop_lost", "exon_loss_variant" ,            
                           "splice_donor_variant", "splice_acceptor_variant"  ); #"5_prime_UTR_premature_start_codon_gain_variant", is considered low impact by snpeff. no need plot
names(arrHighImpactVarTypes) <- gsub('_variant','', arrHighImpactVarTypes);
names(arrHighImpactVarTypes) <- gsub('_',' ', names(arrHighImpactVarTypes) );
#names(arrHighImpactVarTypes)[3] <- "start gain"
 
datVarTypePlot2 <- read.table("../vartype_count_persample.txt", header=T, sep="\t")
datVarTypePlot2 <- rbind(datVarTypePlot2, data.frame(hap=names(arrGeneLossCounts), vartype="ablation", count=arrGeneLossCounts ));

datVarTypePlotHighImpact <- datVarTypePlot2[datVarTypePlot2$vartype %in% arrHighImpactVarTypes, ]

for(i in 1:length(arrHighImpactVarTypes)) {
  datVarTypePlotHighImpact$vartype <- str_replace_all(as.character(datVarTypePlotHighImpact$vartype), as.character(arrHighImpactVarTypes)[i] , names(arrHighImpactVarTypes)[i] )
}

ggplot(datVarTypePlotHighImpact, aes(fill=vartype, y=count, x=hap)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Count', title='High impact variants per haplotype') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(unique(datVarTypePlotHighImpact$vartype)), type = "continuous")) )
ggsave("impacttypes_highimpact_per_hap.pdf", width=6, height=6)

ggplot(datVarTypePlot2, aes(fill=vartype, y=count, x=hap)) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Team', y='Count', title='All variants per haplotype') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_manual('High Impact Variant Type', values=as.character(wes_palette("Zissou1", length(unique(datVarTypePlot2$vartype)), type = "continuous")) )
ggsave("impacttypes_allimpact_per_hap.pdf", width=6, height=6)

#chisq test

mat <- matrix(c(datVarTypePlotHighImpact$count[1:16]), ncol=4, byrow=TRUE )
chisq.test(t(mat) )


mat <- matrix(c(datVarTypePlotHighImpact$count[1:4], datVarTypePlotHighImpact$count[17], 
                datVarTypePlotHighImpact$count[5:8], datVarTypePlotHighImpact$count[18],
                datVarTypePlotHighImpact$count[9:12], datVarTypePlotHighImpact$count[19],
                datVarTypePlotHighImpact$count[13:16], datVarTypePlotHighImpact$count[20]
                ), ncol=5, byrow=TRUE )
chisq.test(t(mat) )
