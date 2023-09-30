#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 32
#SBATCH --mem=120G
CPU=32
genome=/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa


source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate EDTA

ln -sf `realpath $genome` ./genome.fa
RepeatMasker -pa $CPU -s -a -inv -dir ./RepMask -no_is -norna -xsmall -nolow -div 40 -lib ../bahaha.softmasked.fa.mod.EDTA.TElib.fa -cutoff 225 genome.fa

export PERL5LIB=$PERL5LIB:/public/software/conda_envs/EDTA/share/RepeatMasker/

calcDivergenceFromAlign.pl -s RepMask/genome.divsum RepMask/genome.fa.align
./createRepeatLandscape.pl -div RepMask/genome.divsum  -g 691298987  > replandscape.html
