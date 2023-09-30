<?php

$sAssgn = "hap.assign.txt";
$sOutPrefix = "hap.reassigned";
$arrHaps = array("../3ddna/female.hap1.p_ctg.0.review.fasta","../3ddna_hap2/female.hap2.p_ctg.0.review.fasta");
$arrHapFastas = array(fnReadFasta($arrHaps[0]) , fnReadFasta($arrHaps[1]) );
$arrWrittenScfs = array(array(), array());
$sChrPrefix = "bahahaHasm";
$sChrSuffix = "F";
$arrO = array(fopen("$sOutPrefix.0.fa", 'w') , fopen("$sOutPrefix.1.fa", 'w') );

$hA = fopen($sAssgn, 'r');

while(false !== ($sLn = fgets($hA) )) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == 'chr') continue;

        $nChr = intval($arrF[0]);

	foreach(array(0,1) as $nHap) {
		$sHap0Scf =  $arrF[1 + $nHap];

		$sHap0Orient = $arrF[3+$nHap];

		$nHap0To = intval($arrF[13+$nHap]);

		$sScf0Name = "$sChrPrefix"."_".$nChr."_".$sChrSuffix.$nHap0To;

		$sScf0Seq = ($sHap0Orient == '+')? $arrHapFastas[$nHap][$sHap0Scf] : fnRevComp($arrHapFastas[$nHap][$sHap0Scf] );

    echo("Write $sHap0Scf $sHap0Orient to $nHap0To as $sScf0Name\n");
		fwrite($arrO[$nHap0To], ">$sScf0Name\n".wordwrap($sScf0Seq, 100, "\n", true)."\n");

		$arrWrittenScfs[$nHap][$sHap0Scf] = true;
	}
	

}

foreach(array(0,1) as $nHap) {
	foreach(array_keys($arrHapFastas[$nHap]) as $sScf) {
		if (array_key_exists($sScf, $arrWrittenScfs[$nHap])) {
			continue;
		}
    echo("Write unsorted $sScf to $nHap \n");
		fwrite($arrO[$nHap] , ">$sScf\n".wordwrap($arrHapFastas[$nHap][$sScf] , 100, "\n", true)."\n" );
	}
}

function fnReadFasta($sF) {

echo("Loading $sF ...\n");
$arrAsm = array();
$hIn = fopen($sF , 'r');
$sName = '';
$sSeq = '';
do {
        $sLn = fgets($hIn);
        if (false === $sLn || $sLn[0] == '>') {
                if ($sSeq != '') {
                        $arrAsm[$sName] = $sSeq;
                }
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
return $arrAsm;
}


function fnRevComp($sIn) {
        return strrev( strtr($sIn, 'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', 'TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR'));
}


?>
