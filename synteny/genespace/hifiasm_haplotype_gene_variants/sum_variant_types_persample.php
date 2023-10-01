<?php
$sIn = "snpeff_ret.txt"; // count different types of variant counts per sample

$sOut = "vartype_count_persample.txt";

$h = fopen($sIn, 'r');

$arrVartypes = array();
$arrCounts = array(); //key is sample name, subkey is type
$arrSampleIDs = array();
while(false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') {
		continue;
	}

	$arrF = explode("\t", $sLn);
	if ($arrF[0] =='pgChr') {
		$arrSampleIDs = array_slice($arrF, 4);
		$arrSampleIDs = array_slice($arrSampleIDs , 0, count($arrSampleIDs )/2 );
		//print_r($arrSampleIDs );
		//die();
		$arrCounts = array_combine($arrSampleIDs , array_fill(0, count($arrSampleIDs ), array() ) );
		continue;
	}

	$arrF = array_slice($arrF, 4);
	for($i=0;$i<count($arrSampleIDs);$i++) {
		$sSampleID = $arrSampleIDs[$i];
		$arrAnn = fnParseAnnot($arrF[$i]);
		foreach($arrAnn as $sKey => $nVal) {
			if (strpos($sKey, '_variant')!==false || strpos($sKey, 'start_')!==false || strpos($sKey, 'stop_')!==false ) {
				$arrVartypes[$sKey] = true;
				if (!array_key_exists($sSampleID , $arrCounts)) {
					die("$sSampleID not found\n");
				}
				if (!array_key_exists($sKey, $arrCounts[$sSampleID] ) ) {
					$arrCounts[$sSampleID][$sKey] = 0;
				}

				$arrCounts[$sSampleID][$sKey] += intval($nVal);
			}
		}
	}
}

$hO = fopen($sOut, 'w');
fwrite($hO, "hap\tvartype\tcount\n");

foreach( $arrCounts as $sSampleID => $arrVars) {
	foreach($arrVars as $sKey => $nVal) {
		fwrite($hO, "$sSampleID\t$sKey\t$nVal\n");
	}
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
