<?php
$args = getopt("i:o:l:u:");
$sIn = $args['i'];
$sOut = $args['o'];
$nMinRGQ = 30; //min reference hom phred quality
$nMinQUAL = 30;
$nMinGQ = 30;
$nWinSize = 50000; //50kb window
$nStepSize= 100; //250kb step

$nMaxDP = $args['u'];
$nMinDP = $args['l'];

$hIn = popen(" bcftools filter -e ' (GT==\"RR\" && RGQ<$nMinRGQ) || GT==\"mis\" || ( GT!=\"RR\" && (QUAL<$nMinQUAL || GQ<$nMinGQ || INFO/DP>$nMaxDP || INFO/DP< $nMinDP || TYPE!="snp")) ' -O v --threads 2 -S 0 $sIn ", 'r');
$hO = popen("gzip -c > $sOut", 'w'); 

$sScf = '';
$nWinStart = 0;
$nPrevCoord = 0;
$sGenoTrack = ""; //string, 0-hom ref, 1-het, 2-hom alt, 3-missing
while(false !== ($sLn = fgets($hIn) )) {
	$sLn = trim($sLn);
	if ($sLn == '' || $sLn[0] == '#') {
		continue;
	}

	$arrF = explode("\t", $sLn);
	$sLnScf = $arrF[0];
	$nCoord = $arrF[1];
	$sRefBase = strtoupper($arrF[3]);
	$sAltBase = strtoupper($arrF[4]);

	$nStatus = '0'; //default to hom ref.

	if ($sRefBase == 'N' || strpos('ATCG', $sAltBase)==-1 || strlen($sAltBase)!=1 ) {
		$nStatus = '3';
	}

	$arrAnnNames = explode(":", $arrF[8]);
	$arrAnnVals = explode(":", $arrF[9]);
	$arrAnn = array_combine($arrAnnNames, $arrAnnVals);

	if ($arrAnn['GT'] == '0/0' || $arrAnn['GT'] == '0|0') {
		if (array_key_exists('RGQ', $arrAnn) && $arrAnn['RGQ'] < $nMinRGQ) {
			$nStatus = '3';
		} else {
			$nStatus = '0';
		}
	}

	if ($arrAnn['GT'] == '0/1' || $arrAnn['GT'] == '0|1' || $arrAnn['GT'] == '1|0' || $arrAnn['GT'] == '1/0' ) {
		if ($arrF[6]=='PASS') {
			$nStatus = '1';
		} else {
			$nStatus = '3';
		}
	}

	if ($arrAnn['GT'] == '1/1' || $arrAnn['GT'] == '1|1') {

                if ($arrF[6]=='PASS') {
                        $nStatus = '2';
                } else {
                        $nStatus = '3';
                }
	}

        if ($arrAnn['GT'] == './.') {

               $nStatus = '3';
        }

	$sGenoTrack .= $nStatus;

	if ($sScf !== $sLnScf || strlen($sGenoTrack)==($nWinSize) ) {
		//write window
		if (strlen($sGenoTrack) > 100) {
			$sGenoTrack = (strlen($sGenoTrack)==($nWinSize+1))? substr($sGenoTrack, 0, $nWinSize) : $sGenoTrack;
			$nHomRef = substr_count($sGenoTrack, '0');
			$nHet = substr_count($sGenoTrack, '1');
			$nHomAlt = substr_count($sGenoTrack, '2');
			$nMissing = substr_count($sGenoTrack, '3');
			$nEffLen =  $nHomRef + $nHet + $nHomAlt;
			fwrite($hO, "$sScf\t$nWinStart\t".($nPrevCoord+1)."\t".intval(($nWinStart+$nPrevCoord)/2)."\t$nHomRef\t$nHet\t$nHomAlt\t$nMissing\t".($nHet/$nEffLen)."\n" );
		}

		if ($sScf !== $sLnScf) {
			$sScf = $sLnScf;
			$nWinStart = $nCoord;
			$sGenoTrack = $nStatus;
		} else {
			$nWinStart = $nWinStart + $nStepSize;
			$sGenoTrack =  substr($sGenoTrack, $nStepSize ) ;
		}
	}

	$nPrevCoord = $nCoord;
}
?>
