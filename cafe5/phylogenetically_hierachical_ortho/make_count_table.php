<?php
$sIn = "/data/projects/rcui//bahaha_assembly/synteny/genespace/morespp/rundir/orthofinder/Results_Mar18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv";
$sOut = "orthogroup.counts.tsv";

$bHeaderRead = false;

$h = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

while(false !== ($sLn = fgets($h)) ) {
	$sLn = trim($sLn , "\n");
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if (!$bHeaderRead) {
		$arrSpp = array_slice($arrF, 3);
		for($i=0;$i<count($arrSpp);$i++) {
			$arrSpp[$i] = trim($arrSpp[$i]);
		}
		fwrite($hO, "Desc\tOrthogroup\t".implode("\t", $arrSpp )."\n");
		$bHeaderRead = true;
		continue;
	}

	$arrCounts = array();

	for($nCol=3;$nCol<count($arrF);$nCol++) {
		$arrG = explode(',', trim($arrF[$nCol]) ) ;
		$nCount = 0;
		foreach($arrG as $sG) {
			if (trim($sG) !='') {
				$nCount++;
			}
		}
		$arrCounts[] = $nCount;
	}

	fwrite($hO, "(null)\t".trim(str_replace('N0.','', $arrF[0]) )."\t".implode("\t", $arrCounts)."\n");
}

?>
