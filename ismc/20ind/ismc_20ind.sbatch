#!/bin/bash
#SBATCH -p himem
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
sAllTabs=alltabs.tab
gVCF=/public3/group_crf/home/cuirf/bahaha_assembly/wgs/gatk/joined_genotype/joincalled.genotyped.g.vcf.gz
nChr=${SLURM_ARRAY_TASK_ID} 
sOutDir=20inds/${nChr}

sChr=bahahascf_${nChr}

sAllTabs=`realpath $sAllTabs`

mkdir -p $sOutDir
ln -sf `realpath ismc_opt.bpp` $sOutDir/
ln -sf `realpath mapper_opt_100k.bpp` $sOutDir/
ln -sf `realpath mapper_opt_10k.bpp` $sOutDir/

cd $sOutDir

module load bcftools

 bcftools view -M2 -r ${sChr} $gVCF | gzip -c > chr.vcf.gz
 grep -P "$sChr\t" $sAllTabs > in.tab

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate iSMC
ismc params=ismc_opt.bpp
ismc_mapper params=mapper_opt_10k.bpp
ismc_mapper params=mapper_opt_100k.bpp
#rm chr.vcf.gz

exit



