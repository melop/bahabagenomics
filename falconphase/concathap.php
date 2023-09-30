<?php
$nHap = 1;
$sIn = "output/test_data.phased.$nHap.fasta";
$sOut = "concat.hap.$nHap.fasta";

$hIn = fopen($sIn, 'r');
$hOut = fopen($sOut, 'w');
$arrConcat = array();
$sName = '';
$sSeq = '';
do {
	$sLn = fgets($hIn);
	if (false === $sLn || $sLn[0] == '>') {
		fnProcess($sName, $sSeq);
		if (false === $sLn) {
			break;
		}
		$sName = substr(trim($sLn), 1);
		$sSeq = '';
		continue;
	}

	$sSeq .= trim($sLn);
}
while (true);

echo("After concat: ".count( $arrConcat)." contigs \n");

foreach($arrConcat as $sID => $sSeq) {
	fwrite($hOut , ">$sID"."_".$nHap."\n".wordwrap($sSeq, 100, "\n", true)."\n" );
}
function fnProcess($sName, $sSeq) {
	global $arrConcat;
	if ($sName == '') return;
	list($sID) = explode('::', $sName);
	if (!array_key_exists($sID, $arrConcat)) {
		$arrConcat[$sID] = '';
	}

	 $arrConcat[$sID] .= $sSeq;
}

?>
