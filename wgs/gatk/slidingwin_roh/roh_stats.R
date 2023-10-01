args <- commandArgs(T);

dat <- read.table(args[1], header=T);
sSample <- args[2];

#ambstart        homstart        homend  ambend
effGenomeSize <- 648872356;
arrAmbWins <- dat$ambend - dat$ambstart + 1;

nMean <- mean(log10(arrAmbWins));
nSD <-  sd(log10(arrAmbWins));
Froh <- sum(arrAmbWins)/effGenomeSize;

cat(sSample,"\t");
cat(sum(arrAmbWins),"\t" );

cat( nMean ,"\t", nSD, "\t", Froh, "\t");

arrAmbWins <- dat$homend - dat$homstart + 1;

nMean <- mean(log10(arrAmbWins));
nSD <-  sd(log10(arrAmbWins));
Froh <- sum(arrAmbWins)/effGenomeSize;

cat(sum(arrAmbWins),"\t" );

cat( nMean ,"\t", nSD, "\t", Froh, "\n");

