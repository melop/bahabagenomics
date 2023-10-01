ref=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
refname=btp

arrOthers=( larimichthys_crocea collichthys_lucidus nibea_albiflora )
arrFiles=( /data/projects/rcui/bahaha_assembly/annotations/protein_refs/larimichthys_crocea/GCF_000972845.2_L_crocea_2.0_genomic.fna \
/data/projects/rcui/bahaha_assembly/annotations/protein_refs/collichthys_lucidus/GCA_004119915.2_ASM411991v2_genomic.fna \
/data/projects/rcui/bahaha_assembly/annotations/protein_refs/nibea_albiflora/GCA_014281875.1_ASM1428187v1_genomic.fna \
)

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
