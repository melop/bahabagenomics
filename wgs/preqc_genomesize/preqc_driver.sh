for sR1 in ~/../shareddata/bahaha/wgs/trim/BTP*.paired_1.fq.gz; do
        sBase=`basename $sR1`
        sBase=${sBase/.paired_1.fq.gz/}
        sR2=${sR1/.paired_1.fq.gz/.paired_2.fq.gz}
        echo $sBase $sR1 $sR2
        ( bash preqc_BTP.sh $sBase $sR1 $sR2 32 ) > sBase.log 2>&1 
done
wait
