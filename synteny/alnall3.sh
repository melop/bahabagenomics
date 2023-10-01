ref=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
refname=btp

arrOthers=( datnioides_undecimradiatus  lutjanus_erythropterus miichthys_miiuyi sciaenops_ocellatus )
arrFiles=( /data/projects/rcui/bahaha_assembly/annotations/other_refs/Datnioides_undecimradiatus/GCA_008933995.1_BGI_Dund_1.0_genomic.fna \
/data/projects/rcui/bahaha_assembly/annotations/other_refs/Lutjanus_erythropterus/GCA_020091685.1_ASM2009168v1_genomic.fna \
/data/projects/rcui/bahaha_assembly/annotations/other_refs/Miichthys_miiuyi/GCA_001593715.1_ASM159371v1_genomic.fna \
/data/projects/rcui/bahaha_assembly/annotations/other_refs/Sciaenops_ocellatus/GCA_014183145.1_ASM1418314v1_genomic.fna \
)

WD=`pwd`

for i in ${!arrOthers[@]}; do
  cd $WD;
  sp=${arrOthers[$i]}
  file=${arrFiles[$i]}
  dir=$refname.vs.$sp
  mkdir -p $dir
  cd $dir
# ( bash ../aln.sh $ref $file $refname $sp > log.txt 2>&1; bash ../compute_covered.sh  $refname $sp > compute_cov.log 2>&1 ; bash ../paf2pseudogenome.sh $ref $file $refname $sp > pseudogenome.log.txt 2>&1 ) &
  bash ../paf2pseudogenome.sh $ref $file $refname $sp > pseudogenome.log.txt &
done
wait
