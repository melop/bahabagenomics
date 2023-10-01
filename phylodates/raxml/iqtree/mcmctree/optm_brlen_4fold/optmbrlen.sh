ln -s ../../../phylipformat/4fold.phy ./4fold.phy
/data/software/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-AVX2 -t ../../partition.txt.treefile -s 4fold.phy -T 86 -f e  -m GTRGAMMA -n 4foldoptmbrlen 

