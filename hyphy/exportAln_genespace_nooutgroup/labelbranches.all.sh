
arrOutgroupSpp='Sciaenops_ocellatus';
arrUnusedClades='root'; #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren='F'; #c(T , F);

#######################################################################
arrForegroundClades='BTP'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_BTP.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_BTP";

mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

#sOutDIR="Relax_Callopanchax_mrca";
#arrMarkForegroundChildren='F';
#mkdir -p $sOutDIR
#Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

sOutDIR="Codeml_BTP";
arrMarkForegroundChildren='F';
nMarkStyle='codeml';
mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

#sOutDIR="Codeml_Callopanchax_mrca";
#arrMarkForegroundChildren='F';
#mkdir -p $sOutDIR
#Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &
########################################################################

wait
