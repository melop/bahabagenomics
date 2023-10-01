
bedtools annotate -i windows.bed -files repeats.bed repeats.DNA.bed repeats.LTR.bed  cds.norepeats.bed  highimpact.af0-0.2.bed highimpact.af0.2-0.4.bed highimpact.af0.4-0.6.bed \
 highimpact.af0.6-0.8.bed  highimpact.af0.8-0.99.bed  highimpact.af0.99-1.bed  > cov.perwin.bed
