#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/Consurf_slurm-%a.out
#SBATCH -e slurmlog/Consurf_slurm-%a.err
#SBATCH -a 0-19
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load clusterbasics

#SET THIS TO THE CORRECT NUMBER
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate consurf

TOTALPARTS=20

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname


php consurfscores.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

