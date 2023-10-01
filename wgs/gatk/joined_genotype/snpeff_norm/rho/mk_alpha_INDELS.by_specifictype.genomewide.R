setwd("/public3/group_crf/home/cuirf/bahaha_assembly/wgs/gatk/joined_genotype/snpeff_norm/rho");
source("asym.MK.ray.R");
dat <- read.table("consurf_af_sv.w.indelmodifier.rec.bed", header=F, sep="\t")

dat <- dat[ , c(1:3, 6:ncol(dat) ) ]


arrImpact <- c( "frameshift_variant|stop_gained|stop_lost|splice_donor_variant|splice_acceptor_variant|start_lost|exon_loss_variant|gene_fusion|transcript_ablation",
                "upstream_gene_variant|5_prime_UTR_variant", "downstream_gene_variant|3_prime_UTR_variant",  "intron_variant", "intergenic_region"); #"3_prime_UTR_variant", "5_prime_UTR_variant",

arrImpactDesc <- c("High impact variants", "5'-UTR", "3'-UTR", "Intronic", "Intergenic");
lsPlotData <- list();

sOut <- "asymtotoic_mk_indels.byspecifictypes.paired.txt"
arrBreaks <- c(seq(0.0, 0.8, 0.05) , 0.99, 1);
arrAlleleFreq <- arrBreaks[2:(length(arrBreaks)-1)]
AF_low <- 0.3
AF_high <- 0.7
#dat <- read.table("consurf_af_sv.rec.bed", header=F, sep="\t")

dat$AF <- apply((dat[, 10:29]), 1, FUN= function(x) {nTotal = length(x)*2; p <- sum(x)/nTotal; return(p) })
dat$AFBin <- cut(dat$AF, breaks = arrBreaks);
dat <- dat[complete.cases(dat),];

sink(file = sOut)
arrColors <- rev(as.character(wes_palette("Zissou1", length(arrImpact), type = "continuous")))
nCol <- 1;
for(sImpact in arrImpact) {
  #plot asymtotic MK-like :
  dN <- sum(as.integer(grepl(sImpact, dat$V9 ) & dat$AF == 1));
  dS <- sum(as.integer((dat$V10 == 'LOW') & dat$AF == 1));
  pN <- sum(as.integer(grepl(sImpact, dat$V9 ) & dat$AF < 1));
  pS <- sum(as.integer((dat$V10 == 'LOW') & dat$AF < 1));
  
  arrPN <- c();
  arrPS <- c();
  arrAFBins <- as.character(levels(dat$AFBin));
  
  for (sAFBin in arrAFBins[1:(length(arrAFBins)-1)]) {
    arrPN <- c(arrPN, sum(as.integer(grepl(sImpact, dat$V9 ) & dat$AFBin == sAFBin )));
    arrPS <- c(arrPS, sum(as.integer((dat$V10 == 'LOW') & dat$AFBin == sAFBin )));
  }
  

  alpha <- 1 - ( (dS /dN ) * (pN / pS) );
  oConTable <- matrix(c(dN , dS, pN , pS ), nrow = 2);
  oTest <- fisher.test(oConTable);
  cat("\n\n#########################", sImpact, ": #############################\n")
  cat("dN\tdS\tpN\tpS\n");
  cat(dN, "\t", dS, "\t", pN, "\t", pS, "\n");
  cat("alpha global = ", alpha);
  print(oTest)
  
  arrAsymMK <- asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.08, xhigh=0.9, output='table')
  print(arrAsymMK)
  
  n10percIdx <- which(arrAlleleFreq == 0.1);
  alpha_10perc <- 1 - ( (dS * arrPN[n10percIdx] ) / (dN * arrPS[n10percIdx] ) );
  alpha_asymptotic <- arrAsymMK['alpha_asymptotic'];
  
  alphaDiff <- alpha_asymptotic - alpha_10perc;
  
  cat("alpha_10%\talpha_original\talpha_asymptotic\talphaDiff\n");
  cat(alpha_10perc, "\t",  as.numeric(arrAsymMK['alpha_original']) , "\t", as.numeric(alpha_asymptotic), "\t", as.numeric(alphaDiff), "\n");

  lsPlotData[[sImpact]] <- list();
  lsPlotData[[sImpact]]$d0 <- dS;
  lsPlotData[[sImpact]]$d <- dN;
  lsPlotData[[sImpact]]$xlow <- 0.08;
  lsPlotData[[sImpact]]$xhigh <- 0.92;
  lsPlotData[[sImpact]]$df <- data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL);
  lsPlotData[[sImpact]]$true_alpha <- NA;
  lsPlotData[[sImpact]]$arrYLims <- c(-20, 1);
  lsPlotData[[sImpact]]$plotCI <- F;
  lsPlotData[[sImpact]]$col <- arrColors[nCol];
  nCol <- nCol + 1;
}

sink();


# lsPlotData[[1]]$col <- 'red'; #HIGH IMPACT
# lsPlotData[[2]]$col <- 'blue'; #5'utr or upstream
# lsPlotData[[3]]$col <- 'green1'; #3'utr or downstream
# lsPlotData[[4]]$col <- 'yellow1';#intron
# lsPlotData[[5]]$col <- 'violet'; #intergenic
# #lsPlotData[[6]]$col <- 'violet';

pdf(file=paste(sOut, ".pdf",sep="") , width=6, height=4)
#asymptoticMK(d0=dS, d=dN, df=data.frame(f=arrAlleleFreq, p=arrPN, p0=arrPS, row.names=NULL), xlow=0.08, xhigh=0.9, arrYLims=c(-1, 0.5), output='default')

asymptoticMKOverlay(lsPlotData);
#abline(h=alpha_10perc , col='blue');
#abline(h= as.numeric(arrAsymMK['alpha_original'])  , col='green');
#abline(h= 0  , col='black');
legend("bottomright", inset=.05, title="Impact Type", arrImpactDesc,
       lty=1, lwd=2, pch=16, col=arrColors)

dev.off();


