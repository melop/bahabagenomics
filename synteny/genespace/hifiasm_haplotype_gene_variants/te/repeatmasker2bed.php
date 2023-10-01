<?php
$sRepOut = "/public3/group_crf/home/cuirf/bahaha_assembly/annotations/edta/divplot/RepMask/genome.fa.out";
$sOut = "repeats.bed";
$h = fopen($sRepOut, 'r');
$hO = fopen($sOut, 'w');
fgets($h); fgets($h);
while(false!==($sLn = fgets($h) )) {
	$sLn = trim($sLn);

	if ($sLn == '') continue;

	$arrF = preg_split("/\s+/", $sLn);
	$arrO = array_slice($arrF, 4, 3);
	$arrO[] = $arrF[10];
	$arrO[1] = intval($arrO[1] ) - 1;
	fwrite($hO, implode("\t", $arrO)."\n");
}

?>
