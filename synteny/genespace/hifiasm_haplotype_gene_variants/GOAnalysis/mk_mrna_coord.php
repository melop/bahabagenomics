<?php
$sSp = "BTP";
//$sOrthologs = "/data/projects/rcui/bahaha_assembly/annotations/UPhO_morespp/assigngenesymbols/BTP_orthologs_table.ensemblinformed.txt";
$sOrthologs = "/data/projects/rcui/bahaha_assembly/hyphy/exportAln_genespace/genespace.orthogroups.txt";
$sGFF = "/data/projects/rcui/bahaha_assembly/release1.0/annotations_update_symbols/btp.longest.addedUPhO.genesymbol2.spgeneid.gff3";
$sOut = "mRNA.coord.txt";

//OrthoID SpGeneId        GeneSymbol      GeneName        chr     start   end
//Group_360_1     007174  mipb    major intrinsic protein of lens fiber b [Source:ZFIN%3BAcc:ZDB-GENE-050706-86]  chr13   33755782        33756829
//Group_Id        Gene_Symbol     Gene_Name       BTP

$hOrth = fopen($sOrthologs, 'r');
$hGFF = fopen($sGFF, 'r');
$hO = fopen($sOut, 'w');

fwrite($hO, "OrthoID\tspgeneid\tGeneSymbol\tGeneName\tchr\tstart\tend\tmRNAID\n" );
$arrGeneCoord = array();


while(false !== ($sLn = fgets($hGFF) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;

	$arrF = explode("\t", $sLn);
	if (count($arrF)<9) continue;
	//bahahascf_1     maker   mRNA    10095914        10130142        7170.0  -       .       ID=btpa000834-T1;Parent=btpa000834;cds_len=4311;exon_len=4311;exon_insertions=0;exon_deletions=0;sequence=Group_1546_clean_2|larimichthys_crocea|rna-XM_019270470_2;gene_orientation=+;identity=99.30;similarity=99.72;orthoGroup=Group_1546_clean_2;refSp=larimichthys_crocea;refEnsemblProt=rna-XM_019270470_2_hit_1;spgeneid=038404	

	if ($arrF[2] != 'mRNA') continue;

	$sChr = $arrF[0];
	$nStart = $arrF[3];
	$nEnd = $arrF[4];
	preg_match("/spgeneid=([^;]+)/", $arrF[8], $arrM);
	if (count($arrM) != 2 ) continue;

	$sSpGeneID = trim($arrM[1]);
	preg_match("/ID=([^;]+)/", $arrF[8], $arrM);

	$arrGeneCoord[ $sSpGeneID] = array($sChr , $nStart, $nEnd, trim($arrM[1]));
}

echo(count($arrGeneCoord)." spgene id loaded\n" );

$bHeaderRead = false;
$arrHeader = array();

while(false !== ($sLn=fgets($hOrth) )) {

	$sLn = trim($sLn);
	if ($sLn == '') continue;

	$arrF = explode("\t", $sLn);
	if (!$bHeaderRead) {
		$arrHeader = array_flip($arrF);
		if (!array_key_exists($sSp, $arrHeader ) ) {
			die("Error, $sSp not found in ortholog file header!\n");
		}
		$bHeaderRead = true;
		echo("$sSp found on header col ". $arrHeader[$sSp]. "\n");
		continue;
	}

	$sIDString = $arrF[$arrHeader[$sSp]];
	$arrIDs = explode("|", $sIDString);
	$sSpGeneID = $arrIDs[1];
	if (!array_key_exists($sSpGeneID , $arrGeneCoord)) {
		die("spgeneid $sSpGeneID not found in gff!\n");
	}

	fwrite($hO, $arrF[0]."\t".$sSpGeneID ."\t".$arrF[1]."\t".$arrF[2]."\t".implode("\t", $arrGeneCoord[$sSpGeneID] )."\n");
}

?>
