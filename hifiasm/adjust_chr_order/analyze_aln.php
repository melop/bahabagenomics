<?php

$sOut = "hap.assign.txt";
$hO = fopen($sOut, 'w');
$arrAlnMatrix = array();
$nMaxChrNum = 24;
$arrChrs = array();

foreach(array(0,1) as $nHap) {
	$arrAlnMatrix[$nHap] = array();
	foreach(array(0,1) as $nOrigHap) {
		$arrAlnMatrix[$nHap][$nOrigHap] = array();
		$sF = "aln.$nHap.to.hap$nOrigHap.paf";
    $h = fopen($sF, 'r');
    while(false !== ($sLn = fgets($h) ) ) {
      $sLn = trim($sLn);
      if ($sLn == '') continue;
      $arrF = explode("\t", $sLn);
      $sOrigScf = $arrF[0];
      $arrScfName = explode('_', $sOrigScf);
      if (count($arrScfName)!=3) {
        echo("$sOrigScf skip...\n");
        continue;
      }
      list($sScfPrefix, $nScfNum, $sHap) = $arrScfName;
      $nScfNum = intval($nScfNum);
      if ($nScfNum > $nMaxChrNum) {
        continue;
      }

      if (!array_key_exists($nScfNum, $arrAlnMatrix[$nHap][$nOrigHap])) {
        $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum] = array('newscf' => '', 'alnperc' => 0, 'alnorient' => '+' , 'div' => 0);
        $arrChrs[$nScfNum] = true;
      }
            
      $nAlnLen = abs($arrF[3] - $arrF[2]);
      $sAlnOrient = $arrF[4];
      
      $nAlnPerc = 100 * $nAlnLen / $arrF[1];
      if ($nAlnPerc > $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum]['alnperc']) {
        $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum]['newscf'] = $arrF[5];
        $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum]['alnperc'] = $nAlnPerc;
        $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum]['alnorient'] = $sAlnOrient;
        $arrAlnMatrix[$nHap][$nOrigHap][$nScfNum]['div'] = 1-$arrF[9]/$arrF[10];
      }
    }
	}
}

//print_r($arrAlnMatrix);
fwrite($hO, "chr\tscfname0\tscfname1\torient_0\torient_1\talnperc_0_0\talnperc_0_1\talnperc_1_0\talnperc_1_1\tdiv_0_hap0\tdiv_0_hap1\tdiv_1_hap0\tdiv_1_hap1\tassignment_0\tassignment_1\n");
foreach(array_keys($arrChrs) as $nChr ) {
  $arrOut = array();
  $arrOut[] = $nChr;
  $arrOut[] = $arrAlnMatrix[0][0][$nChr]['newscf'];
  $arrOut[] = $arrAlnMatrix[1][0][$nChr]['newscf'];
  $arrOut[] = $arrAlnMatrix[0][0][$nChr]['alnorient'];
  $arrOut[] = $arrAlnMatrix[1][0][$nChr]['alnorient'];
  $arrOut[] = $arrAlnMatrix[0][0][$nChr]['alnperc'];
  $arrOut[] = $arrAlnMatrix[0][1][$nChr]['alnperc'];
  $arrOut[] = $arrAlnMatrix[1][0][$nChr]['alnperc'];
  $arrOut[] = $arrAlnMatrix[1][1][$nChr]['alnperc'];
  $arrOut[] = $arrAlnMatrix[0][0][$nChr]['div'];
  $arrOut[] = $arrAlnMatrix[0][1][$nChr]['div'];
  $arrOut[] = $arrAlnMatrix[1][0][$nChr]['div'];
  $arrOut[] = $arrAlnMatrix[1][1][$nChr]['div'];
  $nAssign0 = ($arrAlnMatrix[0][0][$nChr]['div'] < $arrAlnMatrix[0][1][$nChr]['div'])? 0:1;
  $nAssign1 = ($arrAlnMatrix[1][0][$nChr]['div'] < $arrAlnMatrix[1][1][$nChr]['div'])? 0:1;
  
  if ($nAssign0 == $nAssign1) { //in case initial assignments were the same
    //use secondary standard
    $nMin0 = min($arrAlnMatrix[0][0][$nChr]['div'] , $arrAlnMatrix[0][1][$nChr]['div']); 
    $nMin1 = min($arrAlnMatrix[1][0][$nChr]['div'] , $arrAlnMatrix[1][1][$nChr]['div']); 
    
    if ($nMin0 < $nMin1) {//believe 0
      $nAssign1 = ($nAssign0 == 0)? 1:0;
    } else {
      $nAssign0 = ($nAssign1 == 0)? 1:0;
    }
  }
  $arrOut[] = $nAssign0;
  $arrOut[] = $nAssign1;
  
  
  fwrite($hO, implode("\t", $arrOut)."\n");
  
}

?>
