#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/RelaxNothos_slurm-%a.out
#SBATCH -e slurmlog/RelaxNothos_slurm-%a.err
#SBATCH -a 0-47
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load hyphy/2.5.32
module load clusterbasics

#SET THIS TO THE CORRECT NUMBER
TreeFolder="../exportAln_genespace_nooutgroup/Relax_BTP"
OutDIR=`basename $TreeFolder`;


mkdir -p $OutDIR

php hyphyrelax.php -N ${SLURM_ARRAY_TASK_COUNT} -f ${SLURM_ARRAY_TASK_ID} -T $TreeFolder -o $OutDIR 

exit

