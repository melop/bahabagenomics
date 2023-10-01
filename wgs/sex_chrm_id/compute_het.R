#NOTUSED
setwd("~/bahaha_assembly/wgs/sex_chrm_id/");
library("PopGenome")
GENOME.class <- readVCF("../gatk/BTP_F01.genotyped.g.vcf.gz",numcols=100000, tid="bahahascf_1", from=1, to=34070420, approx=FALSE, out="", parallel=T, gffpath=FALSE, include.unknown = F)

get.sum.data(GENOME.class)
get.neutrality(GENOME.class)

GENOME.class.slide <- sliding.window.transform(GENOME.class,width=100000, jump=50000,type=2,whole.data=T)

GENOME.class.slide <- neutrality.stats(GENOME.class.slide, FAST = T)

ret <- get.neutrality(GENOME.class.slide)

GENOME.class.slide<-count.unknowns(GENOME.class.slide)
GENOME.class@n.unknowns

get.sum.data(GENOME.class.slide)
GENOME.class.slide@n.unknowns

get.biallelic.matrix(GENOME.class)
