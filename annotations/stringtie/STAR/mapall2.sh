#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 32
ulimit -n 99999
module load star/2.7.9a

STAR=STAR


GenomeDir="btpref"
outdir=mapped
tmpdir=tmp
decompressor="zcat"
#GTF=/public/home/shareddata/betta_rnaseq/ref/Betta_splendens.fBetSpl5.2.105.gtf

mkdir -p $outdir
mkdir -p $tmpdir

outdir=`realpath $outdir`
tmpdir=`realpath $tmpdir`

GenomeDir=`realpath $GenomeDir`

for i in ~/../shareddata/bahaha/WHXWZKY-202201122D-01/raw_data_2022-02-24/trim/*.paired_1.fq.gz; do

        sFqFileName=`basename $i`;
	sR2=${i/.paired_1.fq.gz/.paired_2.fq.gz}
        sSampleID=${sFqFileName/.paired_1.fq.gz/};
        sBamFileName=$sSampleID

        echo $sBamFileName;

        rm -rf $tmpdir/${sBamFileName}
        mkdir -p $outdir/${sBamFileName}/
        sDoneFlag=$outdir/${sBamFileName}/all.done;

        if [ -e $sDoneFlag ]; then

                echo $sBamFileName already finished;
        else
                $STAR   --runThreadN $SLURM_CPUS_PER_TASK \
                        --genomeDir "$GenomeDir" \
                        --readFilesIn $i $sR2 \
                        --readFilesCommand "$decompressor" \
                        --outFileNamePrefix $outdir/${sBamFileName}/ \
                        --outSAMmapqUnique 200 \
                        --outSAMtype BAM SortedByCoordinate \
                        --outTmpDir $tmpdir/${sBamFileName}/ \
                        --genomeLoad NoSharedMemory && touch $sDoneFlag 

        fi

done
