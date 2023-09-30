#!/bin/bash
#SBATCH -p himem,hugemem,blade,gpu
#SBATCH -o slurmlog/slurm-%A-%a.out
#SBATCH -e slurmlog/slurm-%A-%a.err
#SBATCH -a 0-63
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#config=simdef2.txt
#config=simdef3.txt
#config=simdef4.txt
#config=simdef5.txt

#config=simdef6.txt
#config=simdef7.txt
#config=simdef8.txt

config=simdef9.txt

module load R


nCounter=0

mkdir -p out

tail -n+2 $config | \
while read Chr TotalGen SeedingPopSize PopSize h h_sd s s_sd RecPerChr TopPercentileForBreeding MeanHomDel SDHomDel MeanDel SDDel GenotypedMarkerPerc Reps; do
    for nThisRep in `seq 1 $Reps`; do
    nModullo=$(( nCounter % SLURM_ARRAY_TASK_COUNT ))
     nCounter=$(( nCounter + 1 ))
     if (( $nModullo == ${SLURM_ARRAY_TASK_ID} )); then
	sCMDLine="Rscript hap_deleterious_offspring_sim_purgeload.cmd.R $nThisRep $Chr $TotalGen $SeedingPopSize $PopSize $h $s $RecPerChr $TopPercentileForBreeding $MeanHomDel $SDHomDel $MeanDel $SDDel $h_sd $s_sd $GenotypedMarkerPerc"
         echo doing job $nCounter $sCMDLine
         $sCMDLine
     fi
    done

done

