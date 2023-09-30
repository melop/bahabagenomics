#!/bin/bash
#SBATCH -p long
#SBATCH -o slurmlog/slurm-%a.out
#SBATCH -e slurmlog/slurm-%a.err
#SBATCH -a 0-107
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1


module load R

arrSeeds=( 1111 )
arSeedingPopSize=( 2 10 100 )
arrPopSize=( 10 100 1000 )
arrH=( "0.5" )
arrS=( "0.003" "0.002" "0.001" "0" )
arrRec=( 2 10 )
arrKeepMinLoadPercentile=( "0.5" "0.25" )

nChr=24
nTotalGen=200;

nCounter=0
for nSeed in "${arrSeeds[@]}"; do
for nSeedingPopSize in "${arSeedingPopSize[@]}"; do
for nPopSize in "${arrPopSize[@]}"; do
for h in "${arrH[@]}"; do
for s in "${arrS[@]}"; do
for nRecBreakpointPerChr in "${arrRec[@]}"; do
for nKeepMinLoadPercentile in "${arrKeepMinLoadPercentile[@]}"; do

	if [ "$s" == "0" ] && [ "$h" != "0" ]; then
    echo skip Rscript hap_deleterious_offspring_sim_plot.cmd.R $nSeed $nChr $nTotalGen $nSeedingPopSize $nPopSize $h $s $nRecBreakpointPerChr
		continue;
	fi
	echo $nCounter
  
	if [ "$nCounter" -eq ${SLURM_ARRAY_TASK_ID} ]; then
	 echo doing hap_deleterious_offspring_sim_purgeload_mutation.cmd.R "$nSeed" "$nChr" "$nTotalGen" "$nSeedingPopSize" "$nPopSize" "$h" "$s" "$nRecBreakpointPerChr" "$nKeepMinLoadPercentile"
   	Rscript hap_deleterious_offspring_sim_purgeload_mutation.cmd.R "$nSeed" "$nChr" "$nTotalGen" "$nSeedingPopSize" "$nPopSize" "$h" "$s" "$nRecBreakpointPerChr" "$nKeepMinLoadPercentile"
    fi
  nCounter=$(( nCounter + 1 ))
done
done
done
done
done
done
done
