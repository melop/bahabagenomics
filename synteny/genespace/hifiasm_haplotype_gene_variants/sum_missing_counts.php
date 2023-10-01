<?php
$sSNPEFF = "/data/software/snpEff/exec/snpeff";
$sGFF2Protein = "/data/software/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl";
$sBEDTOOLS = "bedtools";
$nMaxIntron = 200000; //200kb
$sGenomes = "genomes.txt";
$sSnpEffTablePrefix = "confirmmissing_ret";
$sTmpDir = "confirmmissing_tmp";
$sOutFilePrefix = "counts_deleterious_ret";
$sBLASTDIR = "/data/software/ncbi-blast-2.11.0+/bin/";
$nThreads = 10;
$nAlnPercComplete = 0.9; // if lower than this
$nAlnPercPartial = 0.5; //partial cutoff

echo("reading ...\n");
list( $arrOrthList, $arrSampleIDs) = fnLoadOrthList($sSnpEffTablePrefix) ;
//print_r($arrSampleIDs);
//print_r($arrOrthList);
//die();
//fwrite($hO, "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\t".implode("\t",$arrSampleIDs )."\t".implode("\t",$arrSampleIDs )."\n");

$arrCategories = array('presence', 'HIGH', 'MODERATE', 'LOW', 'MODIFIER');
$arrOutFiles = array();

foreach($arrCategories as $sCat) {
	$arrOutFiles[$sCat] = fopen($sOutFilePrefix.".".$sCat.".txt", 'w');
	fwrite($arrOutFiles[$sCat], "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\t".implode("\t",$arrSampleIDs )."\n");
}

foreach($arrOrthList as $oOrth) {
	$arrSnpEffRet = array();
	$arrMeta = $oOrth['meta'];
	$nOrID = $arrMeta[2];
	list( $sBestSp , $sBestSyn, $sBestProtID) =  explode("|", $oOrth['best']);

	
	foreach($arrCategories as $sCat) {
		
		$arrRetScores = array_fill(0, count($arrSampleIDs), 'NA');
		foreach($arrSampleIDs as $nIdx => $sSample) {
			$oBetterSample = fnChooseBetterGene($oOrth['fields'][$sSample]);
			if ($sCat == "presence") {
				if ($oBetterSample['alnperc']>=$nAlnPercComplete) {
					$arrRetScores[$nIdx] = 1; //treat as complete
				} else if ($oBetterSample['alnperc']>=$nAlnPercPartial) {
					$arrRetScores[$nIdx] = 0.5; //treat as complete
				} else {
					$arrRetScores[$nIdx] = 0; //treat as complete
				}
				continue;
			}

			$arrRetScores[$nIdx] = (array_key_exists($sCat, $oBetterSample)? $oBetterSample[$sCat] : 0 );
		}


		fwrite($arrOutFiles[$sCat], implode("\t", $arrMeta)."\t".implode("\t", $arrRetScores)."\n" );
	}


}

function fnChooseBetterGene($oSample) {
	if (count($oSample)!=2) {
		print_r($oSample);
		die("Error!\n");
	}
	
	foreach(array(0,1) as $nIdx) {
		if (array_key_exists('BEST', $oSample[$nIdx])) {
			$oSample[$nIdx]['alnperc'] = 1;
			return $oSample[$nIdx];
		}
	}

	if (array_key_exists('NA', $oSample[0]) ) {
		$oSample[0]['alnperc'] = 0;
	}

	if (array_key_exists('NA', $oSample[1]) ) {
		$oSample[1]['alnperc'] = 0;
	}

	if ($oSample[0]['alnperc'] > $oSample[1]['alnperc']) {
		return $oSample[0];
	}

	return $oSample[1];
}

function fnLoadOrthList($sPrefix) {
	$hPval = fopen("$sPrefix.txt", 'r');

	$nSampleCount = 0;
	$arrSampleIDs = array();
	$arrCol2ID = array();
	
	$arrMeta = array();
	while(false !== ($sLn = fgets($hPval) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == "pgChr") {
			$nSampleCount = (count($arrF) - 4)/2;
			$arrSampleIDs = array_slice($arrF, 4, $nSampleCount);
			foreach($arrSampleIDs as $nIdx => $sID) {
				$arrCol2ID[4+$nIdx ] = $sID;
				//$arrCol2ID[5+$nSampleCount+$nIdx ] = $sID;
			}
			continue;
		}

		if ($arrF[3] == 'NA') {
			continue;
		}
		$sOrID = $arrF[2];
		$oOrth = array();
		$oOrth['meta'] = array_slice($arrF, 0,4);
		$arrPVals = array();
		foreach($arrCol2ID as $nCol => $sID) {
			$nColSupp = $nCol+$nSampleCount;
			$arrPVals[$sID] = array('NA', 'NA');
			if (array_key_exists($nCol, $arrF)) {
				$arrPVals[$sID][0] = fnParseAnnot($arrF[$nCol]);
				$oBest = isBest($arrF[$nCol]);
				if ($oBest!==false) {
					$oOrth['best'] = "$sID|SYN|".$oBest[0];
					$oOrth['bestP'] = $oBest[1];
				}
			}
			if (array_key_exists($nColSupp, $arrF)) {
				$arrPVals[$sID][1] = fnParseAnnot($arrF[$nColSupp]);
				$oBest = isBest($arrF[$nColSupp]);
				if ($oBest!==false) {
					$oOrth['best'] = "$sID|NS|".$oBest[0];
					$oOrth['bestP'] = $oBest[1];
					
				}
			}

		}
		$oOrth['fields'] = $arrPVals;
		$arrMeta[$sOrID] = $oOrth;
	}


	return array( $arrMeta , $arrSampleIDs);


}

function isBest($s) {
	$arrF = explode('|', $s);
	if ($arrF[0] == 'BEST') {
		if (count($arrF)<2) {
			die($s);
		}
		return explode("=", $arrF[1]);
	}

	return false;
}


function fnParseAnnot($s) {
	$arrF = preg_split('/[|;]/', $s);
	$arrRet = array();
	foreach($arrF as $sF) {
		if (trim($sF) == '') continue;
		$arrFields = explode("=",$sF, 2);
		$sName = "";
		$sVal = 0;
		if (count($arrFields) == 2) {
			list($sName, $sVal) = $arrFields; 
		} else {
			$sName = $arrFields[0];
			$sVal = 1;
		}
		$arrRet[trim($sName, " '\"")] = trim($sVal, " \t'\"");
	}

	return $arrRet;
}

?>