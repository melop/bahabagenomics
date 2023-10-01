ref=/data/projects/rcui/bahaha_assembly/annotations/repeatmask/scf.softmasked.fa
refname=btp

arrOthers=( cheilinus_undulatus gasterosteus_aculeatus micropterus_salmoides sander_lucioperca anabas_testudineus astyanax_mexicanus \
betta_splendens scleropages_formosus channa_argus gadus_morhua oryzias_latipes neogobius_melanostomus poecilia_formosa takifugu_rubripes tetraodon_nigroviridis \
oreochromis_niloticus danio_rerio \
)
arrFiles=( /data/projects/rcui//bahaha_assembly/annotations/other_refs/Cheilinus_undulatus/GCF_018320785.1_ASM1832078v1_genomic.fna \
/data/projects/rcui//bahaha_assembly/annotations/other_refs/Gasterosteus_aculeatus/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna \
/data/projects/rcui//bahaha_assembly/annotations/other_refs/Micropterus_salmoides/GCF_014851395.1_ASM1485139v1_genomic.fna \
/data/projects/rcui//bahaha_assembly/annotations/other_refs/Sander_lucioperca/GCF_008315115.2_SLUC_FBN_1.2_genomic.fna \
/data/projects/rcui//macropodus_compare/ensembl/anabas/Anabas_testudineus.fAnaTes1.2.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/astyanax/Astyanax_mexicanus.Astyanax_mexicanus-2.0.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/betta/Betta_splendens.fBetSpl5.2.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/bonytongue/Scleropages_formosus.fSclFor1.1.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/channa_argus/GCA_004786185.1_ASM478618v1_genomic.fna \
/data/projects/rcui/macropodus_compare/ensembl/gadus/Gadus_morhua.gadMor3.0.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/medaka/Oryzias_latipes.ASM223467v1.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/neogobius/Neogobius_melanostomus.RGoby_Basel_V2.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/poecilia/Poecilia_formosa.PoeFor_5.1.2.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/takifugu/Takifugu_rubripes.fTakRub1.2.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/tetraodon/Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/tilapia/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna_sm.toplevel.fa \
/data/projects/rcui/macropodus_compare/ensembl/zebrafish/Danio_rerio.GRCz11.dna_sm.toplevel.fa \
)

WD=`pwd`

for i in ${!arrOthers[@]}; do
  cd $WD;
  sp=${arrOthers[$i]}
  file=${arrFiles[$i]}
  dir=$refname.vs.$sp
  mkdir -p $dir
  cd $dir
 ( bash ../aln.sh $ref $file $refname $sp > log.txt 2>&1; bash ../compute_covered.sh  $refname $sp > compute_cov.log 2>&1 ; bash ../paf2pseudogenome.sh $ref $file $refname $sp > pseudogenome.log.txt 2>&1 ) &
#  bash ../paf2pseudogenome.sh $ref $file $refname $sp > pseudogenome.log.txt &
done
wait
