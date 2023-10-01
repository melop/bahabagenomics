setwd("/data/projects/rcui/bahaha_assembly/wgs/gatk/joined_genotype/popkin_relatedness");
library(BEDMatrix)
library(popkin)

X <- BEDMatrix("converted.bed")
dim(X)
kinship <- popkin(X)
plot_popkin(
  kinship, 
  # shared bottom and left margin value, to make space for labels
  mar = 1
)

hist(kinship);
write.table(kinship, file = "kinship.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
#inbr(kinship)

write.table(inbr_diag(kinship), file = "kinship_diaginbrdcoef.tsv", quote = F, row.names = T, col.names = T, sep = "\t")


pairwise_fst <- pwfst(kinship)

leg_title <- expression(paste('Pairwise ', F[ST]))
# NOTE no need for inbr_diag() here!
plot_popkin(
  pairwise_fst,
  labs_even = TRUE,
  labs_line = 1,
  labs_cex = 0.7,
  leg_title = leg_title,
  mar = c(2, 0.2)
)
