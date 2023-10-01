ref=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
refname=btp

arrOthers=( larimichthys_crocea_jimeiU )
arrFiles=( /data/projects/rcui/bahaha_assembly/annotations/protein_refs/larimichthys_crocea_jimeiU/GCA_004352675.2_ASM435267v2_genomic.fna )
WD=`pwd`

for i in ${!arrOthers[@]}; do
  cd $WD;
  sp=${arrOthers[$i]}
  file=${arrFiles[$i]}
  dir=$refname.vs.$sp
  mkdir -p $dir
  cd $dir
#  bash ../aln.sh $ref $file $refname $sp > log.txt 2>&1 &
#  bash ../compute_covered.sh  $refname $sp > compute_cov.log 2>&1 &
bash ../paf2pseudogenome.sh $ref $file $refname $sp > log.txt 2>&1 &
done
wait
