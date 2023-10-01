<?php
//map pgID to zebrafish and human Ensembl IDs, for GO analysis
$sIDMap = "../hifiasm_bahaha_haplotypes_zebrafish_human_forGO/synorthos.txt"; //bahaha transcript IDs to zebrafish and human
$arrEnsemblSpp = array("BahahaTaipingensis"); //the ensembl species to use in $sIDMap, if multiple genes, choose the first one
$sThisMap = "orth_length_outliers.ids.txt";
$sOut = "pgid2BTPRef.txt";
//first load mapping between bahaha ids and zebrafish/human ids
$h = fopen($sIDMap, 'r');

$arrSampleCols = array();
$arrRefCols = array();
$nSampleRefCount = 0;
$arrMapTranscriptID2RefEnsemblID = array();

while(false !== ($sLn = fgets($h)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);

	if ($arrF[0] == "pgChr") {
		$arrSampleIDs = array_slice($arrF, 3);
		$arrSampleIDs = array_slice($arrSampleIDs, 0, count($arrSampleIDs)/2 );
		$nSampleRefCount = count($arrSampleIDs);
		$arrSampleIDFlip = array_flip($arrSampleIDs);
		foreach($arrEnsemblSpp as $sEnsemblSp) {
			if (!array_key_exists($sEnsemblSp, $arrSampleIDFlip)) {
				die("Error: ensembl species $sEnsemblSp not found in $sIDMap \n");
			}

			$arrRefCols[$sEnsemblSp] = $arrSampleIDFlip[$sEnsemblSp]+3;
			unset($arrSampleIDFlip[$sEnsemblSp]);
		}

		foreach($arrSampleIDFlip as $sSampleID => $nCol) {
			$arrSampleCols[$sSampleID] = $nCol + 3;
			$arrMapTranscriptID2RefEnsemblID[$sSampleID] = array();
		}

		//print_r($arrSampleCols);
		//print_r($arrRefCols);
		//print_r($nSampleRefCount);
		continue;
	}

	$arrRefSpIDs = array();

	foreach($arrRefCols as $sRefID => $nRefCol) {
		$nRefNSCol = $nRefCol + $nSampleRefCount;
		if ($arrF[$nRefCol] == 'NA' && $arrF[$nRefNSCol] == 'NA') {
			$arrRefSpIDs[$sRefID] = 'NA';
			continue;
		}

		if ($arrF[$nRefCol] != 'NA') {
			$arrRefIDs = explode(";", $arrF[$nRefCol]);
			$arrRefSpIDs[$sRefID] = $arrRefIDs[0];
			continue;
		}

		if ($arrF[$nRefNSCol] != 'NA') {
			$arrRefIDs = explode(";", $arrF[$nRefNSCol]);
			$arrRefSpIDs[$sRefID] = $arrRefIDs[0];
			continue;
		}

	}

	foreach($arrSampleCols as $sSampleID => $nSampleCol) {
		$nSampleNSCol = $nSampleCol + $nSampleRefCount;
		if ($arrF[$nSampleCol] == 'NA' && $arrF[$nSampleNSCol] == 'NA') {
			continue;
		}
		$arrTransIDs = array();
		if ($arrF[$nSampleCol] != 'NA') {
			$arrTransIDs = array_merge($arrTransIDs , explode(';',$arrF[$nSampleCol] ));
		}

		$arrNSTransIDs = array();
		if ($arrF[$nSampleNSCol] != 'NA') {
			$arrNSTransIDs = array_merge($arrNSTransIDs , explode(';',$arrF[$nSampleNSCol] ));
		}

		foreach($arrTransIDs as $sTransID) {
				$arrMapTranscriptID2RefEnsemblID[$sSampleID][$sTransID] = $arrRefSpIDs;
		}

		foreach($arrNSTransIDs as $sTransID) { //non-syn ids are secondary choice
			if (!array_key_exists($sTransID,$arrMapTranscriptID2RefEnsemblID[$sSampleID]) ) {
				$arrMapTranscriptID2RefEnsemblID[$sSampleID][$sTransID] = $arrRefSpIDs;
			}
		}

	}
}

//print_r($arrMapTranscriptID2RefEnsemblID);
$h2 = fopen($sThisMap, 'r');
$nSampleCount = 0;
$arrSampleCols = array();
$hO = fopen($sOut, 'w');
fwrite($hO, "pgid\t".implode("\t", array_keys($arrRefCols))."\n");

while(false !== ($sLn = fgets($h2)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == "pgChr") {
		$arrSampleIDs = array_slice($arrF, 5);
		$arrSampleIDs = array_slice($arrSampleIDs, 0, count($arrSampleIDs)/2 );
		$nSampleCount = count($arrSampleIDs );
		$arrSampleIDFlip = array_flip($arrSampleIDs);
		foreach($arrSampleIDFlip as $sSampleID => $nCol) {
			$arrSampleCols[$sSampleID] = $nCol + 5;
		}

		//print_r($arrSampleCols);
		//print_r($nSampleCount);
		//die();
		continue;
	}

	$arrRet = array();
	foreach($arrRefCols as $sRefID => $nTmp) {
		$arrRefGeneIDCounts = array();
		foreach($arrSampleCols as $sSampleID => $nSampleCol) {
			$nSampleNSCol = $nSampleCol + $nSampleCount;

			$sTransID = $arrF[$nSampleCol];
			$sNSTransID = $arrF[$nSampleNSCol];

			if ($sTransID != 'NA') {

				$sMappedTo = "NA";
				if (array_key_exists($sTransID, $arrMapTranscriptID2RefEnsemblID[$sSampleID])) {
					$sMappedTo = $arrMapTranscriptID2RefEnsemblID[$sSampleID][$sTransID][$sRefID];
				}
				if ($sMappedTo != 'NA') {
					
					if (!array_key_exists($sMappedTo, $arrRefGeneIDCounts)) {
						$arrRefGeneIDCounts[$sMappedTo] = 0;
					}

					$arrRefGeneIDCounts[$sMappedTo]  += 10;
					continue;
				}
			}

			if ($sNSTransID != 'NA') {

				$sMappedTo = "NA";
				if (array_key_exists($sNSTransID, $arrMapTranscriptID2RefEnsemblID[$sSampleID])) {					
					$sMappedTo = $arrMapTranscriptID2RefEnsemblID[$sSampleID][$sNSTransID][$sRefID];
				}

				if ($sMappedTo != 'NA') {
					
					if (!array_key_exists($sMappedTo, $arrRefGeneIDCounts)) {
						$arrRefGeneIDCounts[$sMappedTo] = 0;
					}

					$arrRefGeneIDCounts[$sMappedTo]  += 1;
					continue;
				}
			}

		}

		//print_r($arrRefGeneIDCounts);
		if (count($arrRefGeneIDCounts)==0) {
			$arrRet[$sRefID ] = 'NA' ;
			continue; //no correspondence found
		}
		$nMaxScore = max($arrRefGeneIDCounts);
		$arrRefGeneIDCountsFlip = array_flip($arrRefGeneIDCounts);
		$sMappedID = $arrRefGeneIDCountsFlip[$nMaxScore];
		$arrRet[$sRefID ] = $sMappedID ;
	}

	fwrite($hO, $arrF[2]."\t".implode("\t",$arrRet)."\n");

}
?>
