<?php
$sSnpEffVcf = "ann.merged.vcf.gz";
$sCloseOutgroup = "Miichthys_miiuyi";
$sFarOutgroup = "larimichthys_crocea_jimeiU";
$sOut = "consurf.SV.polarized.out.txt";

$hO = fopen($sOut, 'w');

//print_r($arrConsurf);

$hVCF = popen("zcat -f $sSnpEffVcf", 'r');

$nCloseOutgroupCol = 0;
$nFarOutgroupCol = 0;
$arrHaps = array();
$arrCovMask = array();
while(false !== ($sLn = fgets($hVCF) )) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == "#CHROM") {
		//read in header line
		$arrRev = array_flip($arrF);
		if (!array_key_exists($sCloseOutgroup, $arrRev)) {
			die("Outgroup sp. $sCloseOutgroup not found in vcf header!\n");
		} else {
			$nCloseOutgroupCol = $arrRev[$sCloseOutgroup];
			unset($arrRev[$sCloseOutgroup]);
		}

		if (!array_key_exists($sFarOutgroup, $arrRev)) {
			die("Outgroup sp. $sFarOutgroup not found in vcf header!\n");
		} else {
			$nFarOutgroupCol =  $arrRev[$sFarOutgroup];
			unset($arrRev[$sFarOutgroup]);
		}

		$arrHaps = array_slice($arrRev, 9);
		//$arrCovMask = fnLoadMask(array_keys($arrHaps));
		$arrMaskSpp = array_keys($arrHaps );
		$arrMaskSpp[] = $sCloseOutgroup;
		$arrMaskSpp[] = $sFarOutgroup;
		$arrCovMask = fnLoadMask($arrMaskSpp );
		//print_r($arrHaps);
		//echo("$nFarOutgroupCol $nCloseOutgroupCol\n");
		//die();
		continue;

	} 
	if ($sLn[0] == '#') {
		continue;
	}

	// if they differ, then the ancestral state cannot be ascertained
	if ( (!fnInRegion($sCloseOutgroup, $arrF[0] , $arrF[1])) || (!fnInRegion($sFarOutgroup, $arrF[0] , $arrF[1])) ) { //outgroups must have coverage at this site
		continue;
	}
	if ($arrF[$nCloseOutgroupCol] != $arrF[$nFarOutgroupCol]) {
		continue;
	}


	$arrAnn = fnParseSnpEffAnn($arrF[7]);
	if (count($arrAnn) == 0) {
		continue;
	}

	
	//found overlapped consurf score and variants
	//check to make sure that both outgroups are the same genotype
	$sAncestralState = $arrF[$nCloseOutgroupCol];
	$sAncAllele = $arrF[3];
	$sDerAllele = $arrF[4];
	if ($sAncestralState == "1/1") {
		$sAncAllele = $arrF[4];
		$sDerAllele = $arrF[3];
	} else if ($sAncestralState != "0/0" ) {
		continue;
	}

	$arrPolarizedSamples = array(); 

	foreach($arrHaps as $sSample => $nCol) {
		if (!fnInRegion($sSample, $arrF[0] , $arrF[1])) {
			$arrPolarizedSamples[] = 'NA';
		} else if ($arrF[$nCol] == $sAncestralState) {
			$arrPolarizedSamples[] = 0;
		} else {
			$arrPolarizedSamples[] = 1;
		}
	}


	fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$sAncAllele."\t".$sDerAllele."\t".implode("\t", $arrAnn[0]) ."\t".implode("\t", $arrPolarizedSamples)."\n");

}

function fnLoadMask($arrNames) {
	$arrRet = array();
	foreach($arrNames as $sName) {
		$sF = "pairwise.simp.ref.$sName.filteredhomologchr.txt";
		if (!file_exists($sF)) {
			$sF = "pairwise.simp.ref.".$sName.".txt";
		}
		echo("Loading mask $sF ...\n");
		if (!file_exists($sF)) {
			die("Mask $sF not found.\n");
		}

		$h = fopen($sF, 'r');
		$arrRet[$sName] = array();
		while(false !== ($sLn = fgets($h))) {
			$sLn = trim($sLn);
			$arrF = explode("\t", $sLn);
			$sScf = $arrF[4];
			$nStart = $arrF[5];
			$nEnd = $arrF[6];
			if (!array_key_exists($sScf, $arrRet[$sName])) {
				$arrRet[$sName][$sScf] = array();
			}

			$arrRet[$sName][$sScf][] = array($nStart, $nEnd );
		}

		foreach(array_keys($arrRet[$sName]) as $sScf) {
			usort($arrRet[$sName][$sScf], 'fnSortCoord');
		}
	}


	return $arrRet;
}

function fnSortCoord($arrCoord1, $arrCoord2) {
	$a = min($arrCoord1);
	$b = min($arrCoord2);
	if ($a == $b) {
        return 0;
    }
    return ($a < $b) ? -1 : 1;
}

function fnInRegion($sSample, $sScf, $nPos) {
	global $arrCovMask;
	if (!array_key_exists($sScf, $arrCovMask[$sSample])) {
		return false;
	}

	foreach($arrCovMask[$sSample][$sScf] as $arrCoords) {
		if ($nPos > $arrCoords[0] && $nPos <= $arrCoords[1]) {
			return true;
		}

	}

	return false;
}

function fnParseSnpEffAnn($s) {
	$arrF1 = explode(";", $s);
	$arrMissenseVars = array();
	foreach($arrF1 as $sField) {
		list($sKey, $sVal) = explode('=', $sField);
		if ($sKey != 'ANN') {
			continue;
		}

		$arrAnns = explode(',', $sVal);
		foreach($arrAnns as $sAnnStr) {
			$arrAnnInfo = explode('|', $sAnnStr);
			if ($arrAnnInfo[1]=='missense_variant' ||  $arrAnnInfo[2]=='MODIFIER' ) {
				continue;
			}
			$sRNAID = $arrAnnInfo[6];
			$sProteinVar = $arrAnnInfo[10];
			$arrMissenseVars[] = array( $sRNAID , $arrAnnInfo[1], $arrAnnInfo[2], $sProteinVar );
		}
	}

	return $arrMissenseVars;
}
function fnLoadConsurf($s) {
	$h = popen("zcat -f $s", "r");
	$arrRet = array();
	while(false !== ($sLn = fgets($h) )) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		$sCoord = $arrF[0].":".$arrF[1];
		if (!array_key_exists($sCoord, $arrRet)) {
			$arrRet[$sCoord] = array();
		} 

		$arrRet[$sCoord][$arrF[3]] = $arrF[2];
	}

	return $arrRet;
}
?>
