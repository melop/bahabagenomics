outstem=joincalled
genome=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa 
gatk=/data/software/gatk-4.2.1.0/gatk
CPU=36 #12
SCRATCHDIR=`pwd`/tmp
vcfs="../BTP_???.g.vcf.gz ../BTP_Q0.g.vcf.gz ../BTP_Q1.g.vcf.gz"

ulimit -n 99999 #change file handler limit!

mkdir -p $SCRATCHDIR



#STEP 6: combine gvcf
STEP=06.combinegvcf.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    ls $vcfs > gvcfs.list
    $gatk --java-options "-Xmx180G -XX:ParallelGCThreads=4" CombineGVCFs -R $genome --variant gvcfs.list -O $outstem.combined_gvcf.vcf > $STEP.log 2>&1 \
    && touch  $STEP.done

    bgzip -@ $CPU $outstem.combined_gvcf.vcf >> $STEP.log 2>&1 && tabix -p vcf $outstem.combined_gvcf.vcf.gz >> $STEP.log 2>&1

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#STEP 7: genotype gvcf
STEP=07.genotypegvcf.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    $gatk --java-options '-Xmx120g' GenotypeGVCFs -O $outstem.genotyped.g.vcf -R $genome -V $outstem.combined_gvcf.vcf.gz  -all-sites true > $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi

#STEP 8: compress
STEP=08.compress.$outstem

if [ -f $STEP.done ]; then
    echo "$STEP done, skip...";
else
    bgzip -@ $CPU $outstem.genotyped.g.vcf > $STEP.log 2>&1 && tabix -p vcf $outstem.genotyped.g.vcf.gz >> $STEP.log 2>&1 \
    && touch  $STEP.done

    if [  "$?" -ne  "0" ]; then
      exit 1;
    fi
fi
