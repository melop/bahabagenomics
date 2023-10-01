<?php
$sRunDir = "rundir";

$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Bahaha_taipingensis_F0";
$sVersion = "1.0";
$sLongestIsoformProt = "/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate1.8.15/hifiasm_hap_F0/btpf0.longest_isoform.prot.fa";

//$sSp = "Bahaha_taipingensis_F1";
//$sVersion = "1.0";
//$sLongestIsoformProt = "/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate1.8.15/hifiasm_hap_F1/btpf1.longest_isoform.prot.fa";


//$sSp = "Bahaha_taipingensis_M0";
//$sVersion = "1.0";
//$sLongestIsoformProt = "/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate1.8.15/hifiasm_hap_M0/btpm0.longest_isoform.prot.fa";

//$sSp = "Bahaha_taipingensis_M1";
//$sVersion = "1.0";
//$sLongestIsoformProt = "/public3/group_crf/home/cuirf/bahaha_assembly/annotations/funannotate1.8.15/hifiasm_hap_M1/btpm1.longest_isoform.prot.fa";



$sWD = "$sRunDir/rawGenomes/$sSp/$sVersion/annotation";
exec("mkdir -p $sWD");

$hGFF = false;
if ($sScfPrefix !='') {
	$hGFF = popen("sort -k1,1n -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
} else {
	$hGFF = popen("sort -k1,1 -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
}

$hFas = popen("gzip -c > $sWD/$sSp.pep.fa.gz", 'w');

$hIn = popen("zcat -f $sLongestIsoformProt", 'r');

while(false!== ($sLn = fgets($hIn) )) {
	$sLn = trim($sLn);
	if ($sLn == '') {
		continue;
	}

	if ($sLn[0] == '>') {
		$arrF = explode(' ', substr($sLn, 1));
		$sSeqName = $arrF[0];
		$sSeqName = preg_replace('/[:|]/', '_', $sSeqName);
		$sPos = $arrF[count($arrF)-1];
		preg_match('/([^:]+):([^-]+)-([^(]+)\\((\\S)\\)/', $sPos, $arrM);
		if (count($arrM) != 5) {
			die("Failed to parse header: $sLn\n");
		}

		list($sChr, $nStart, $nEnd, $sOrient) = array_slice($arrM, 1);
		$sChr = str_replace($sScfPrefix, '', $sChr);
		fwrite($hGFF, "$sChr\tfunannotate\tgene\t$nStart\t$nEnd\t\t$sOrient\t\tlocus=$sSeqName\n");
		fwrite($hFas, ">$sSeqName\n");
	} else {
		$sLn = str_replace('*', '', $sLn);
		fwrite($hFas, $sLn."\n");
	}

	
}

?>
