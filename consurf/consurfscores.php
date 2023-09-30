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

$sRefSp = "REF";

$sRefGenome = "/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/genome/bahaha.softmasked.fa";
$sGTFFile = "/public3/group_crf/home/cuirf/bahaha_assembly/release1.0/annotations_orthofinder_symbols/btp.longest.addedUPhO.genesymbol2.spgeneid.gff3";

/*fnParseConsurfRet("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/consurf/consurf_run/659392274/out/PLP_consurf.grades", "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/consurf/consurf_run/659392274/pos.map.txt");
die();
*/
while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-o':
            $sOutDIR  = trim(array_shift($argv));
            break;
        case '-S':
            $sGeneSymbolMap  = trim(array_shift($argv));
            break;
        case '-N':
            $nTotalPart  = trim(array_shift($argv));
            break;
        case '-f':
            $nThisPart  = trim(array_shift($argv));
            break;

    }
}



exec("mkdir -p $sOutDIR");
$sOutFile = "$sOutDIR/ret_part$nThisPart.of.$nTotalPart.txt";
$arrComputedResults = fnLoadCurrentRet($sOutFile ); //load previous results, and redo items that failed or not done yet.

$hOut = fopen($sOutFile , 'w'); //append to it


//list($arrGeneSymbolMap,$arrSpGeneId2GroupId) = fnLoadGeneSymbolMap($sGeneSymbolMap);
list($arrGeneSymbolMap,$arrSpGeneId2GroupId) = array(array(), array());


//============================================================================
$hRefGenome = fopen($sRefGenome, "r");
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




$sLn = "";

$sSeqName = "";
$sSeq = "";
echo("Loading genome fasta...\n");
do {
	$sLn = fgets($hRefGenome);

	if ($sLn==false) {
		//fnProcessRecord($sSeqName, $sSeq);
		$arrRefGenome[$sSeqName] = $sSeq;
		break;
	}

	$sLn = trim($sLn);
	if (substr($sLn , 0,1) == ">") {
		$arrRefGenome[$sSeqName] = $sSeq;
		$sSeqName = substr($sLn , 1);
		$sSeq = "";
		
	} else {
		$sSeq .= $sLn;
	}

} while (true);




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
	$sGroupId = "NA";

	if (!array_key_exists($sSpGeneId, $arrSpGeneId2GroupId)) {
		echo("warning: spgeneid $sSpGeneId not found in ortholog table\n");
		//continue;
	} else {
		$sGroupId = $arrSpGeneId2GroupId[$sSpGeneId];
	}

	echo("spgeneid = $sSpGeneId, groupid = $sGroupId \n");

	$bOnReverse = ($arrRNA["strand"]=="-");
	$sScfld = $oGene["scf"];
	$arrExclusionList = array();// no need to exclude fnBuildExclusionList($arrRNA,  array('3primeuncertain', '5primeuncertain', 'hugeinsertionlist' )); // exclude the AA positions potentially caused by mis-annotation.

	if (count($arrExclusionList) > 0) {
		echo("Exclude following AA positions due to unreliable gene annotation $arrRNAIDs[0] :  ");
			foreach($arrExclusionList as $key => $value) {
			  echo("$key-$value ");
			}
		echo("\n");
	}

	//sort:
	if ($bOnReverse) {
		krsort($arrRNA['CDSs']);
	} else {
		ksort($arrRNA['CDSs']);
	}

	$nPrevStart = -1;
	$nPrevEnd = -1;
	$sPrevLeftOverBases = ""; //this is the left over bases at the end of an exon ,to be joined with the next exon.
	$nAccuAALen = 0; //this is the accumulative AA position on the protein.
	$sFilteredCDSSeq = "";
	$arrGenomicPos = array();

	foreach($arrRNA['CDSs'] as $nExonStart => $nExonEnd ) { 


		//now go over each exon, note that if orientation is + , then exon starts should monotonically increase, otherwise decrease
		if ($nPrevStart != -1) {
			if ( ($bOnReverse && $nPrevStart <$nExonStart) ||  (!$bOnReverse && $nPrevStart > $nExonStart)) {
				print_r($arrRNA);
				die("Error, in the above gene the exon order is wrong!");
			}
		}



		if (!array_key_exists($sScfld , $arrRefGenome) ) {
			die("CONTIG $sScfld not found in genome fasta.\n");
		}
		$sSeqExon = substr($arrRefGenome[$sScfld],$nExonStart-1, $nExonEnd -  $nExonStart + 1);
		if ($bOnReverse ) {
			$sSeqExon = fnReverseComplement($sSeqExon);
		}

		$nBpAppendedFromPrev = strlen($sPrevLeftOverBases); // appended this many bases at the begining.
		$sSeqExon = $sPrevLeftOverBases . $sSeqExon;
		$nLeftOverBases = strlen($sSeqExon) % 3;
		$sLeftOverBases = substr($sSeqExon , strlen($sSeqExon) - $nLeftOverBases); #bug fixed on 21.7.2016
		$nInFrameBases = strlen($sSeqExon) - $nLeftOverBases;
		for($i=0;$i<$nInFrameBases;$i+=3) {
			$nAccuAALen++;

			if (fnIsExcluded($nAccuAALen , $arrExclusionList)) {
				//echo($nAccuAALen." AA excluded. $arrRNAIDs[0]\n");
				//print_r($arrExclusionList);
				//die();
				continue;
			}

			$sCodon = substr( $sSeqExon, $i, 3);

			//loop over the 3 codon bases:
			for ($nCodonOffset=0;$nCodonOffset < 3; $nCodonOffset++) {
				$nPos = fnGetPos($i + $nCodonOffset, strlen($sPrevLeftOverBases), $bOnReverse, $nPrevStart, $nPrevEnd , $nExonStart, $nExonEnd ); //returns 1-based coordinates
				
				$arrGenomicPos[] = $nPos;
			}

			$sFilteredCDSSeq .= $sCodon;



		}

		$sPrevLeftOverBases = $sLeftOverBases;
		$nPrevStart = $nExonStart;
		$nPrevEnd = $nExonEnd;
	}

	if (strlen($sFilteredCDSSeq) != count($arrGenomicPos) ) {
		echo("Length mismatch between cleaned CDS and positions\n");
		continue;
	}

	//print_r($sFilteredCDSSeq);
	//print_r($arrGenomicPos);
	fnProcessAln($sSpGeneId, $sGroupId , array('REF' => $sFilteredCDSSeq), $arrGenomicPos);
	//die();

} 

//============================================================================

function fnProcessAln($sGeneName, $sGroupId , $arrFullCodingSeq, $arrGenomicPos) {
	global $sOutDIR, $nThisPart, $nTotalPart, $nGeneCount,$hOut,$arrComputedResults,$arrGeneSymbolMap, $sRefSp, $sConsurfTmp;

	$nGeneCount++;
	if ($nGeneCount % $nTotalPart != $nThisPart) return;

	if (array_key_exists($sGeneName, $arrComputedResults)) {
		echo("$sGeneName done, skip\n");
		fwrite($hOut, implode("\t", $arrComputedResults[$sGeneName] )."\n");
		return;
	}

	$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);

	$arrTaxa = array_keys($arrFullCodingSeq);

	if (!array_key_exists($sRefSp , $arrFullCodingSeq)) {
		echo("Ref species $sRefSp not found in alignment , skip:\n");
		echo(implode("\t",$arrTaxa)."\n");
		return;
	}



	echo("Gene : $sGeneName\n");
	$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.


		if (!$arrAlnRet["stopcodon"] ) { //there is no stop codon, do codeml

			echo(" ... Doing Consurf... \n");
			$sWD = realpath(fnGetConsurfTmp());
			$sConsurfOut = "$sWD/out";
			$sProbeProteinRaw = "$sWD/probeAARaw.fa";
			$sProbeProtein = "$sWD/probeAA.fa";
			$sPosMap = "$sWD/pos.map.txt";
			$hProbeProteinRaw = fopen($sProbeProteinRaw , 'w');
			$hProbeProtein = fopen($sProbeProtein , 'w');
			$hPosMap = fopen($sPosMap , 'w');			
			$oAAStrippedMissing = fnStripMissing($arrAlnRet['proteins'][$sRefSp]);

			//write position map:
			foreach($oAAStrippedMissing['map'] as $nNewCoord => $nOldCoord) {
				fwrite($hPosMap , "$nNewCoord\t$nOldCoord\n");
			}

			fwrite($hProbeProteinRaw , ">$sRefSp\n".$arrAlnRet['proteins'][$sRefSp]);
			fwrite($hProbeProtein , ">$sRefSp\n".$oAAStrippedMissing['AA']);
			$sCmd = "mkdir -p $sConsurfOut; cd $sWD ; consurf -Out_Dir $sConsurfOut -m -Seq_File probeAA.fa -MaxHomol 150  -Algorithm LikelihoodML --Iterat 1";
			echo($sCmd."\n");
			exec($sCmd );
			$sConsurfGrades = "$sConsurfOut/$sRefSp"."_consurf.grades";
			if (!file_exists($sConsurfGrades) ) {
				echo("Consurf failed!\n");
                                fwrite($hOut, "$sGeneName\tFailed\t$sWD\n");
				return;
			}

			$arrConsurfRet = fnParseConsurfRet($sConsurfGrades, $sPosMap);
			$arrConsurfRetGenomicPos = fnMatch2Genomic($arrConsurfRet, $arrGenomicPos);
                        fwrite($hOut, "$sGeneName\t$sGroupId\tSuccess\t");
			foreach($arrConsurfRet  as $nPos => $nScore) { 
				fwrite($hOut, "$nPos:$nScore;");
			}

                        fwrite($hOut, "\t");

			foreach($arrConsurfRetGenomicPos  as $nPos => $nScore) { 
				fwrite($hOut, "$nPos:$nScore;");
			}

                        fwrite($hOut, "\t$sConsurfOut\n");

			
		}
		else {
			echo("$sGeneName\tStopCodonFoundIn:".implode(",",$arrAlnRet["stoptaxa"]).PHP_EOL);
		}


	
}

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

