<?php
$hIn = popen("cat out/ret_part* | grep 'Success'", 'r');
$hOut = popen("gzip -c > consurf_scores_flat.txt.gz", 'w');
$sMap = "spgeneid_map_coord.txt";

$arrMap = fnParseMap($sMap);
while( false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	$arrFF = explode(";", $arrF[4]);
	foreach($arrFF as $sVal) {
		$arrFFF = explode(':', $sVal);
		if (count($arrFFF)==2) {
			fwrite($hOut , $arrMap[$arrF[0]]."\t$arrFFF[0]\t$arrFFF[1]\t$arrF[0]\n");
		}
	}
}

function fnParseMap($sMap) {
	$h = fopen($sMap, 'r');
	$arrRet = array();	
	while( false !== ($sLn = fgets($h)) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		if (count($arrF) <3) continue;
		$arrRet[$arrF[0]] = $arrF[2];
	}

	return $arrRet;
	
}

?>
