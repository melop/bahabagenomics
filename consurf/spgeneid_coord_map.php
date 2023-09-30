<?php
//This script takes the exported alignments and compute the consurf score per amino acid position
//using homologous proteins recruited from the uniprot database.
//PLP is taken as the reference sequence because it is an outgroup
require_once(dirname(__FILE__) . "/lib.php");

$sOutDIR = "out";
$sConsurfTmp = "consurf_run";
//$sGeneSymbolMap = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";

$nThisPart = 0;
$nTotalPart = 1;


//$sRefGenome = "/beegfs/group_dv/data/human_GRCh37.p13/GRCh37.p13_chr_rename.fa";
$sGTFFile = "/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/annotations_orthofinder_symbols/btp.longest.addedUPhO.genesymbol2.spgeneid.gff3";

$sOutFile = "spgeneid_map_coord.txt";

$hOut = fopen($sOutFile , 'w'); //append to it


//list($arrGeneSymbolMap,$arrSpGeneId2GroupId) = fnLoadGeneSymbolMap($sGeneSymbolMap);


//============================================================================
$arrRefGenome = array(); //save ref genome into the RAM



	$sGTF = $sGTFFile;//ENSXMAG00000001050";//ENSXMAG00000001400";//ENSXMAG00000000500//ENSXMAG00000000995//ENSXMAG00000001272

	if (!file_exists($sGTF)) {
		echo("Warning: $sGTF not found\n");
		die();
	}
	else {
		echo("Loading GTF: $sGTF ...\n");
		//die();
	}

	$oGTF = new MappedCDSGFF();//XiphoGTF();
	$oGTF->LoadGFF($sGTF, '', 'CDS'); //use "exon" features




//now go over the GTF file:
$arrGenes = $oGTF->GetGenes('maker');

foreach($arrGenes as $oGene) {
	$arrRNAIDs = array_keys($oGene['mRNAs']);
	if (count($arrRNAIDs) ==0 ) {
		continue;
	}
	$arrRNA = $oGene['mRNAs'][$arrRNAIDs[0]]; // only use the first mRNA

	//print_r($arrRNA);

	if ( (!array_key_exists('ID' , $arrRNA['annot'])) || trim($arrRNA['annot']['ID'])=='' ) {
		continue;
	}

	$sSpGeneId = $arrRNA['annot']["ID"];
	$sGroupId = $arrRNA['annot']["ID"];

	echo("spgeneid = $sSpGeneId, groupid = $sGroupId \n");

	$bOnReverse = ($arrRNA["strand"]=="-");
	$sScfld = $oGene["scf"];
	fwrite($hOut , "$sSpGeneId\t$sGroupId\t$sScfld\t".$arrRNA["strand"]."\t".$arrRNA["start"]."\t".$arrRNA["end"]."\n");
} 

//============================================================================


function fnMatch2Genomic($arrConsurfRet, $arrGenomicPos) {
	$arrRet = array();
	foreach($arrConsurfRet  as $nAAPos => $nScore) { //$nAAPos is 0 based, the keys of $arrGenomicPos is also 0 based nucleotides
		$nPos1 = $nAAPos * 3 + 0;
		$nPos2 = $nAAPos * 3 + 1;
		$nPos3 = $nAAPos * 3 + 2;
		$arrRet[$arrGenomicPos[$nPos1]] = $nScore;
		$arrRet[$arrGenomicPos[$nPos2]] = $nScore;
		$arrRet[$arrGenomicPos[$nPos3]] = $nScore; //this may be synonymous, but still return sites. 
	}

	return $arrRet;
}

function fnGetConsurfTmp() {
	global $sConsurfTmp;
	if (!file_exists($sConsurfTmp)) {
		exec("mkdir -p $sConsurfTmp");
	}
	$sRet = "";
	do {
		$sRet = $sConsurfTmp . "/" .rand();
	} while(file_exists($sRet ) );
	exec("mkdir -p $sRet");

	return $sRet;

}

function fnStripMissing($s) {
	$nLen = strlen($s);
	$sStripProtein = "";
	$arrPosMap = array(); //key is new position in the stripped protein, value is the original position on the gapped sequence

	for($i=0;$i<$nLen; $i++) {
		if ($s[$i] == 'X' || $s[$i] == '-' || $s[$i] == '?' ) {
			continue;
		}		

		$sStripProtein .= $s[$i];
		$arrPosMap[] = $i;
	}

	return array('map' => $arrPosMap , 'AA' => $sStripProtein);
}

function fnParseConsurfRet($sConsurfGrades, $sPosMap) {
	$hPosMap = fopen($sPosMap , 'r');
	$hConsurfGrades = fopen($sConsurfGrades , 'r');

	$arrPosMap = array();
	while( false !== ($sLn = fgets($hPosMap) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$arrPosMap[$arrF[0]] = $arrF[1];
	}

	//print_r($arrPosMap);

	$arrScores = array();
	while( false !== ($sLn = fgets($hConsurfGrades) ) ) {

		if ($sLn == '') continue;
		if ($sLn[0] == "\t") continue;
		if ($sLn[0] == "-") continue;
		if ($sLn[0] == "*") break;

		$arrF = preg_split("/\s+/", trim($sLn) );
		if ( count( $arrF) < 8) continue;
		if ($arrF[0] == 'POS') continue;
		if (strpos($arrF[3] , '*') !== false ) continue;
		$nNewPos = $arrF[0] - 1;
		if (!array_key_exists($nNewPos , $arrPosMap) ) {
			echo("Aln pos $nNewPos not found in positional map\n");
			continue;
		}

		$arrScores[$arrPosMap[$nNewPos]] = $arrF[3];
	}

	//print_r($arrScores);
	//die();
	return $arrScores;
}

function fnLoadCurrentRet($sOutFile ) {

	$arrRet = array();

	if (!file_exists($sOutFile) ) {
		echo("Previous output $sOutFile not found, new run\n ");
		return $arrRet;
	}

	$hPrevOut = fopen($sOutFile , 'r');

	while( false !== ($sLn = fgets($hPrevOut) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$arrRet[$arrF[0]] = $arrF;
	}
	fclose($hPrevOut);

//	$hOut = fopen($sOutFile , 'w');

	foreach($arrRet as $sGroup => $arrLn) {
		if ($arrLn[1] == 'Success'  ) {
//			fwrite($hOut , implode("\t" , $arrLn) . "\n"); //only write back good ones.
		} else {
			exec("rm -R $arrLn[2]");
			echo("To be retried: $sGroup \n");
		}
	}

//	fclose($hOut);

	return $arrRet;
}

function fnLoadGeneSymbolMap($sOrthologList, $sSp="NOR") {

	if (!file_exists($sOrthologList) ) {
		die("Gene symbol definition not found $sOrthologList\n");
	}

	$hOrthologList = fopen($sOrthologList , "r");
	$arrOrthologMeta = array();
	$arrSpGeneId2GroupId = array();

	$nLn = -1;
	$nSpCol = 0;
	echo("Parsing ortholog definitions...\n");
	while( ($sLn=fgets($hOrthologList))!==false ) {
		$sLn = trim($sLn);
		if ($sLn == "") {
			continue;
		}
		$nLn++;
		if ($nLn==0) {
			$arrF = explode("\t", $sLn);
			$arrFFlip = array_flip($arrF);
			if (array_key_exists($sSp, $arrFFlip)) {
				$nSpCol = $arrFFlip[$sSp];
			}

			continue; // skip header
		}
	
		$arrFields = explode("\t", $sLn);
		$arrOrthologMeta[$arrFields[0]] = array_slice($arrFields , 1, 2);
		$arrSpGeneId2GroupId[$arrFields[$nSpCol]] = $arrFields[0];
	}

	echo("Loaded ". count($arrOrthologMeta) ." ortholog definitions\n" );

	return array($arrOrthologMeta , $arrSpGeneId2GroupId);
}

function str_lreplace($search, $replace, $subject)
{
    $pos = strrpos($subject, $search);

    if($pos !== false)
    {
        $subject = substr_replace($subject, $replace, $pos, strlen($search));
    }

    return $subject;
}

function fnGetPos($i , $nPrevLeftoverBases , $bOnReverse, $nPrevStart, $nPrevEnd , $nExonStart, $nExonEnd ) {
	//check if this is the bases from the previous exon:

	if ( $i  < $nPrevLeftoverBases ) { // return index based on previous exon
		if ($bOnReverse) {
			return $nPrevStart+ ($nPrevLeftoverBases -1) - $i;
		} else {
			return $nPrevEnd - ($nPrevLeftoverBases-1) + $i;
		}
	} else {
		$i = $i - $nPrevLeftoverBases ;
		if ($bOnReverse) {
			return $nExonEnd - $i;
		} else {
			return $nExonStart + $i;
		}
	}
}

function fnBuildExclusionList($arrRNA,  $arrExcludeTypes ) {
	$nAALen = $arrRNA['CDSLen'] / 3;
	$arrRet = array();
	//print_r($arrExcludeTypes);
	//print_r($arrRNA);
	foreach($arrExcludeTypes as $sType) {
		if ( (!array_key_exists($sType , $arrRNA['annot'])) || trim($arrRNA['annot'][$sType])=='' ) {
			continue;
		}

		if ($sType == '5primeuncertain') {
			if ($arrRNA['annot'][$sType] > 0) {
				$arrRet[1] = $arrRNA['annot'][$sType];
			}
			continue;
		}
		if ($sType == '3primeuncertain') {
			if ($arrRNA['annot'][$sType] > 0) {
				$arrRet[$nAALen-$arrRNA['annot'][$sType]+1] = $nAALen;// bug fix, 2018-07-9. this bug may cause masking of whole gene.
			}
			continue;
		}

		$arrList = explode(',', $arrRNA['annot'][$sType]);
		if (count($arrList) > 0 ) {
			foreach($arrList as $sRange) {
				
				 $arrSplitRange = explode('-' , $sRange);
				if (count( $arrSplitRange)!=2) continue;
				list($nExcStart, $nExcEnd)  = $arrSplitRange;
				$arrRet[$nExcStart] = $nExcEnd;
			}
		}
		
	}
    ksort($arrRet);	
    return $arrRet;
}

function fnIsExcluded($nAccuAALen , &$arrExclusionList) {
	if (count($arrExclusionList) == 0 ) return false;

	foreach($arrExclusionList as $nExcStart => $nExcEnd) {
		if ($nExcStart <= $nAccuAALen && $nAccuAALen <= $nExcEnd) {
			return true;
		}
	}

	return false;
	
}



?>

