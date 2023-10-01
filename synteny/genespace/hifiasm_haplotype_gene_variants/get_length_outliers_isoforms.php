<?php
ini_set('memory_limit', '50000M');

//perform pairwise blast on outgroup peptides, get a distribution of aligned lengths
//then blast each ingroup protein against the outgroup proteins, get the mean length, and compare this mean to the outgroup alignment distribution, and mark any ingroup peptide that looks abnormal
$sTmp = "diamond_tmp";
$nThreads = 24;
$sDIAMOND = "diamond";
$sOut = "orth_length_outliers.txt";
$sOut2 = "orth_length_outliers.ids.txt";
$nMinOutgroup = 3;
$nMinStdev = 0.01; //used when computed stddev is less than this ; in the unit of fraction. for example 0.98 +- 0.01

$sInList = "../hifiasm_bahaha_haplotypes/synorthos.txt";
$sProteinDir = "../hifiasm_bahaha_haplotypes/rundir/peptide/"; //longest isoform proteins

$sIngroupList = "ingroups.txt";
$sOutgroupList = "outgroup.txt";

$arrIngroups = fnLoadIngroupList($sIngroupList);
$arrOutgroups = fnLoadList($sOutgroupList);

$arrProteins = fnLoadProteins( $arrIngroups , $arrOutgroups, $sProteinDir );
$arrIsoformMaps = fnLoadIsoformIDMap($arrIngroups);
//print_r($arrIngroups);
//print_r($arrOutgroups);

$descriptorspec = array(
   0 => array("pipe", "r"),  // 标准输入，子进程从此管道中读取数据
   1 => array("pipe", "w"),  // 标准输出，子进程向此管道中写入数据
   2 => array("pipe", "w") // 标准错误，写入到一个文件
);

$process = proc_open('Rscript pnorm.R', $descriptorspec, $PNormPipes);


$hIn = fopen($sInList, 'r');
list($hO, $hO2, $arrProcessedRecords, $bAppendMode) = fnCheckPrevResults($sOut, $sOut2);
//print_r($arrProcessedRecords);
//die();

$arrIngroupCols = array();
$arrOutgroupCols = array();
while(false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == 'pgChr') {
		//process first line
		$nTotalSpp = (count($arrF)-3)/2;
		$arrFFlip = array_flip( array_slice($arrF , 0, $nTotalSpp+3) );

		$arrIngroupCols = array_intersect_key($arrFFlip, array_flip(array_keys($arrIngroups)) );
		$arrOutgroupCols = array_intersect_key($arrFFlip, array_flip($arrOutgroups) );

		//print_r($arrIngroupCols);
		//print_r($arrOutgroupCols);

		//die();
                if ($bAppendMode) {
                        //append to previous result files, no need to write header
                        echo("Previous output files exist, resuming from break point ...\n");
			continue;
                }

		fwrite($hO, "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\tOutgroupAlnPercStdev\t".implode("\t",array_keys($arrIngroups) )."\t".implode("\t",array_keys($arrIngroups) )."\n");
		fwrite($hO2, "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\tOutgroupAlnPercStdev\t".implode("\t",array_keys($arrIngroups) )."\t".implode("\t",array_keys($arrIngroups) )."\n");
		continue;

	}

	$sCheckIDStr = $arrF[0]."_".$arrF[1]."_".$arrF[2];
	if (array_key_exists($sCheckIDStr, $arrProcessedRecords)) {
		echo("skip $sCheckIDStr \n");
		continue;
	}
	//process:
	//check the outgroup sequences
	echo("Processing $sCheckIDStr \n");
	$arrOutgroupSeqs = array();
	foreach($arrOutgroupCols as $sSp => $nCol) {
		if ($arrF[$nCol] == "NA") {
			continue;
		}
		$sProtID = $arrF[$nCol];
		if (!array_key_exists($sProtID, $arrProteins[$sSp])) {
			die("Error 1! protein $sProtID not found in species $sSp\n");
		}
		$arrOutgroupSeqs["$sSp|$sProtID"] = $arrProteins[$sSp][$sProtID];
	}

	if (count($arrOutgroupSeqs) <$nMinOutgroup ) {
		fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$arrF[2]."\tNA\tNA\t".implode("\t", array_fill(0, count($arrIngroups) , 1) )."\n");
		continue;
	}
	//now perform alignment:
	exec("mkdir -p $sTmp");
	$sFAS = "$sTmp/in.fas";
	$hFAS = fopen($sFAS, 'w');
	foreach($arrOutgroupSeqs as $sID => $sSeq) {
		fwrite($hFAS, ">$sID\n$sSeq\n");
	}
	fclose($hFAS);
	$bExecStatus = exec("cd $sTmp; rm outgroup*; $sDIAMOND makedb --in in.fas -d outgroup; $sDIAMOND blastp -p $nThreads -d outgroup -q in.fas -o outgroup.csv --outfmt 6 qseqid sseqid qlen slen positive mismatch gapopen gaps evalue > screen.log 2>&1");
	//die();
	if ($bExecStatus === false ) {
		die("Diamond command failed!\n");
	} 
	$arrBlastRet = fnReadBlastRet("$sTmp/outgroup.csv");
	$arrSeqIDs = array_keys($arrOutgroupSeqs);
	$arrAlnPerc = array();
	foreach($arrSeqIDs as $sID1) {
		foreach($arrSeqIDs as $sID2) {
			if ($sID1 == $sID2) {
				continue;
			}
			$sPairString = "$sID1 vs. $sID2";
			if (!array_key_exists($sPairString,$arrBlastRet )) {
				echo("Warning: $sPairString blast hit not found in result!\n");
				continue;
			}
			$arrAlnPerc[] = $arrBlastRet[$sPairString];
		}
		
	}

	if (count($arrAlnPerc)==0) {
		fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$arrF[2]."\tNA\tNA\t".implode("\t", array_fill(0, count($arrIngroups) , 1) )."\n");
		continue;		
	}

	if (count($arrAlnPerc)==1) {
		fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$arrF[2]."\t".$arrAlnPerc[0]."\tNA\t".implode("\t", array_fill(0, count($arrIngroups) , 1) )."\n");
		continue;		
	}

	//print_r($arrAlnPerc);

	//die();

	list($nMeanAlnLen, $nStddevAlnLen) = mean_stddev($arrAlnPerc);
	echo("mean $nMeanAlnLen, stddev $nStddevAlnLen\n");
	//die();
	//next check each ingroup sequence against outgroups:
	$arrIngroupSeqs = array();
	$arrIngroupRet = array();
	$arrIngroupProtID = array();
	foreach($arrIngroupCols as $sSp => $nCol) {
		$arrIngroupRet[$sSp] = 1;
		if ($arrF[$nCol] == "NA") {
			$arrIngroupRet[$sSp] = -1; //set to -1 if missing
			$arrIngroupProtID[$sSp] = 'NA' ;
			continue;
		}
		$sProtLongestID = $arrF[$nCol];

		//TODO: add isoform selection code
		if (!array_key_exists($sProtLongestID, $arrIsoformMaps[$sSp])) {
			die("Error 2! Representative protein $sProtLongestID not found in isoform map, please check!\n");
		}

		$arrPValRet = array();
		foreach($arrIsoformMaps[$sSp][$sProtLongestID] as $sProtID) {
			echo("Checking isoform for $sProtLongestID : $sProtID ...\n");
			if (!array_key_exists($sSp, $arrProteins)) {
				die("Error 3! species $sSp not found in loaded proteins\n");
			}
			if (!array_key_exists($sProtID, $arrProteins[$sSp])) {
				die("Error 4! protein $sProtID not found in species $sSp\n");
			}
			file_put_contents("$sTmp/q.fas", ">$sSp|$sProtID\n".$arrProteins[$sSp][$sProtID]."\n");
			$bExecStatus = exec("cd $sTmp; $sDIAMOND blastp -p $nThreads -d outgroup -q q.fas -o q.csv --outfmt 6 qseqid sseqid qlen slen positive mismatch gapopen gaps evalue > screen.log 2>&1");
			if ($bExecStatus === false ) {
				die("Diamond command failed!\n");
			}
			$arrBlastRet = fnReadBlastRet("$sTmp/q.csv");
				//print_r($arrBlastRet );
				//die();
			$nMeanAlnPerc = 0;
                        $nPVal = 0;
                        if (count($arrBlastRet)>0) {
                        	$nMeanAlnPerc = (float)array_sum($arrBlastRet)/((float)count($arrBlastRet));
                                $nPVal = PNorm($nMeanAlnPerc,$nMeanAlnLen, $nStddevAlnLen);
                        } else {
				echo("Warning: diamond produced no blast results.\n");
			}
                        echo("mean = $nMeanAlnPerc, p = $nPVal\n");

			$arrPValRet[intval($nPVal*100000)] = array('p' => $nPVal, 'alnperc' => $nMeanAlnPerc, 'prot_id' => $sProtID);
		}
		$arrBestIsoform = $arrPValRet[max(array_keys($arrPValRet))];
		echo("Selected best isoform: ".$arrBestIsoform['prot_id']."\t".$arrBestIsoform['p']."\t".$arrBestIsoform['alnperc']."\n");
		$arrIngroupRet[$sSp] = $arrBestIsoform['p'];
		$arrIngroupProtID[$sSp] = $arrBestIsoform['prot_id'] ;
	}

	//check other orthologs:
	$arrIngroupSuppRet = array();
	$arrIngroupSuppID = array();

	foreach($arrIngroupCols as $sSp => $nCol) {
		$arrIngroupSuppRet[$sSp] = -1;
		$nNSOrCol = $nCol + $nTotalSpp;
		if ($arrF[$nNSOrCol] == "NA") {
			$arrIngroupSuppRet[$sSp] = -1; //set to -1 if missing
			$arrIngroupSuppID[$sSp] = 'NA';
			continue;
		}	

		$arrProtIDs = explode(";", $arrF[$nNSOrCol]);
		$nMaxAlnPercP = 0;
		$sMaxAlnID = 'NA';
		foreach($arrProtIDs as $sProtLongestID) {
			
		        if (!array_key_exists($sProtLongestID, $arrIsoformMaps[$sSp])) {
        	                die("Error 5! Representative protein $sProtLongestID not found in isoform map, please check!\n");
	                }
			foreach($arrIsoformMaps[$sSp][$sProtLongestID] as $sProtID) {
				if ($sMaxAlnID == 'NA') {
					$sMaxAlnID = $sProtID;
				}
				if (!array_key_exists($sProtID, $arrProteins[$sSp])) {
					die("Error 6! protein $sProtID not found in species $sSp\n");
				}
				file_put_contents("$sTmp/q2.fas", ">$sSp|$sProtID\n".$arrProteins[$sSp][$sProtID]."\n");
				$bExecStatus = exec("cd $sTmp; $sDIAMOND blastp -p $nThreads -d outgroup -q q2.fas -o q2.csv --outfmt 6 qseqid sseqid qlen slen positive mismatch gapopen gaps evalue > screen.log 2>&1");
				if ($bExecStatus === false ) {
					die("Diamond command failed!\n");
				}
				$arrBlastRet = fnReadBlastRet("$sTmp/q2.csv");
					//print_r($arrBlastRet );
					//die();
				if (count($arrBlastRet) == 0) {
					$arrIngroupSuppRet[$sSp] = -1; //set to -1 if missing
					//die("blastnotfound\n");
					continue;
				}
				$nMeanAlnPerc = 0;
				$nPVal = 0;
				if (count($arrBlastRet)>0) {
					$nMeanAlnPerc = (float)array_sum($arrBlastRet)/((float)count($arrBlastRet));
					$nPVal = PNorm($nMeanAlnPerc,$nMeanAlnLen, $nStddevAlnLen);
				} else {
        	                        echo("Warning: diamond produced no blast results.\n");
	                        }

				echo("mean = $nMeanAlnPerc, p = $nPVal\n");
				$nMaxAlnPercP = max($nMaxAlnPercP, $nPVal );
				if ($nPVal == $nMaxAlnPercP)  {
					$sMaxAlnID = $sProtID;
				}
			}
		}

		$arrIngroupSuppRet[$sSp] = $nMaxAlnPercP;
		$arrIngroupSuppID[$sSp] = $sMaxAlnID;
	}
	//print_r($arrIngroupRet);
	//die();
	fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$arrF[2]."\t$nMeanAlnLen\t$nStddevAlnLen\t".implode("\t", $arrIngroupRet)."\t".implode("\t", $arrIngroupSuppRet)."\n");
	fwrite($hO2, $arrF[0]."\t".$arrF[1]."\t".$arrF[2]."\t$nMeanAlnLen\t$nStddevAlnLen\t".implode("\t", $arrIngroupProtID)."\t".implode("\t", $arrIngroupSuppID)."\n");

}



fclose($PNormPipes[0]);
fclose($PNormPipes[1]);
fclose($PNormPipes[2]);
proc_close($process);

function fnReadBlastRet($s) {
	$h = fopen($s, 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($h)) ) {
			$sLn = trim($sLn);
			if ($sLn == '') continue;
			$arrF = explode("\t", $sLn);
			$sProt1 = $arrF[0];
			$sProt2 = $arrF[1];
			if ($sProt1 == $sProt2) {
				continue;
			}
			$nProt1Len = $arrF[2];
			$nProt2Len = $arrF[3];
			$nAlnLen = $arrF[4];

			$nPercAln1 = $nAlnLen / $nProt1Len;
			$sPairString1 = "$sProt1 vs. $sProt2";
			$nPercAln1 = (array_key_exists($sPairString1, $arrRet)? max($arrRet[$sPairString1], $nPercAln1) : $nPercAln1);

			$nPercAln2 = $nAlnLen / $nProt2Len;
			$sPairString2 = "$sProt2 vs. $sProt1";
			$nPercAln2 = (array_key_exists($sPairString2, $arrRet)? max($arrRet[$sPairString2], $nPercAln2) : $nPercAln2);

			$arrRet[$sPairString1] = $nPercAln1;
			$arrRet[$sPairString2] = $nPercAln2;
	}

	return $arrRet;
}
function PNorm($x, $mean, $stddev) {
	global $PNormPipes;
	fwrite($PNormPipes[0], "$x\t$mean\t$stddev\n");
	return floatval(fgets($PNormPipes[1]) );
}

function fnLoadIngroupList($sF) {
        $arrRet = array();
        $h = fopen($sF, 'r');
        while(false !== ($sLn = fgets($h)) ) {
                $sLn = trim($sLn);
                if ($sLn != '') {
			$arrF = explode("\t", $sLn);
			if (count($arrF)!=3) {
				echo("Warning: expecting two fields in $sF, $sLn\n");
				continue;
			}
                        $arrRet[$arrF[0]] = array('isoformlist' => $arrF[1] , 'proteinfile' => $arrF[2]);
                }
        }

        return $arrRet;

}

function fnLoadList($sF) {
	$arrRet = array();
	if (!file_exists($sF)) {
		die("$sF not found! \n");
	}
	$h = fopen($sF, 'r');
	while(false !== ($sLn = fgets($h)) ) {
		$sLn = trim($sLn);
		if ($sLn != '') {
			$arrRet[] = $sLn;
		}
	}

	return $arrRet;
}

function fnLoadProteins($arrIngroup, $arrOutgroup, $sProteinDir) {
	$arrRet = array();
	foreach($arrOutgroup as $sSp) {
		$sFile = "$sProteinDir/$sSp.fa";
		if (!file_exists($sFile)) {
			die("Error 7! protein file $sFile not found\n");
		}

		echo("Loading $sFile ...\n");
		$arrRet[$sSp] = readFasta($sFile);
	}

        foreach($arrIngroup as $sSp => $arrInfo) {
                $sFile = $arrInfo['proteinfile'];
                if (!file_exists($sFile)) {
                        die("Error 8! protein file $sFile not found\n");
                }

                echo("Loading $sFile ...\n");
                $arrRet[$sSp] = readFasta($sFile);
        }


	return $arrRet;
}

function fnLoadIsoformIDMap($arrIngroup) {

        foreach($arrIngroup as $sSp => $arrInfo) {
                $sFile = $arrInfo['isoformlist'];
                if (!file_exists($sFile)) {
                        die("Error 9! isoform list file $sFile not found\n");
                }

                echo("Loading $sFile ...\n");
                $arrRet[$sSp] = readIsoformList($sFile);
        }

	return $arrRet;
}

function readIsoformList($s) {
        $h = fopen($s, 'r');
        $arrRet = array();
	echo("Loading isoform map $s ...\n");
	while(false!==($sLn = fgets($h))) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 2) {
			continue;
		}
		$arrIsoforms = array_slice($arrF, 1);
		foreach($arrIsoforms as $sIsoform) {
			$arrRet[$sIsoform] = $arrIsoforms;
		}
	}

	return $arrRet;
}

function readFasta($s) {
        $h = fopen($s, 'r');
        $arrRet = array();
        $sName = '';
        $sSeq = '';
        do {
                $sLn = fgets($h);
                if (false===$sLn) {
                        if ($sSeq!='') {
                                $arrRet[$sName] = trim($sSeq ,'*');

                        }
                        break;
                }

                $sLn = trim($sLn);
                if ($sLn == '') continue;
                if ($sLn[0]=='>') {
                        if ($sSeq!='') {
                                $arrRet[$sName] = trim($sSeq , '*');

                        }
                        $sName = substr($sLn, 1);
			list($sName) = preg_split("/\s+/", $sName);
                        $sSeq = '';
                        continue;
                }

                $sSeq .= $sLn;

        } while(true);

        return $arrRet;
}

function mean_stddev($arr)
    {
        $num_of_elements = count($arr);
          
        $variance = 0.0;
          
                // calculating mean using array_sum() method
        $average = array_sum($arr)/$num_of_elements;
          
        foreach($arr as $i)
        {
            // sum of squares of differences between 
                        // all numbers and means.
            $variance += pow(($i - $average), 2);
        }
          
        return array($average , (float)sqrt($variance/$num_of_elements) );
    }

//list($hO, $hO2, $arrProcessedRecords, $bAppendMode) = fnCheckPrevResults($sOut, $sOut2);
function fnCheckPrevResults($sOut, $sOut2) {
	if ((!file_exists($sOut)) || (!file_exists($sOut2)) ) {
		return array(fopen($sOut,'w'), fopen($sOut2, 'w'), array(), false );
	}

	$h1 = fopen($sOut, 'r');

	$arrDone = array();
	while(false !== ($sLn = fgets($h1) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);

		$sIDStr = implode("_", array_slice($arrF, 0,3));
		$arrDone[$sIDStr] = true;
	}
	fclose($h1);
	return array(fopen($sOut,'a'), fopen($sOut2, 'a'), $arrDone, true );
}

?>
