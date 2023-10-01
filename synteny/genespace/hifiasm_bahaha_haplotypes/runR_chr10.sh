source /opt/miniconda3/etc/profile.d/conda.sh

export PATH=$PATH:/data/software/MCScanX:/opt/miniconda3/bin/
#conda activate orthofinder
Rscript run_chr10.R > run10.log 2>&1
