<?php
ini_set('memory_limit', '50000M');

$sSNPEFF = "/public/software/snpEff/exec/snpeff";
$sGFF2Protein = "/public/software/conda_envs/funannotate1.8.15/opt/pasa-2.5.3/misc_utilities/gff3_file_to_proteins.pl";
$sBEDTOOLS = "/public/software/conda_envs/funannotate1.8.15/bin/bedtools";
$nMaxIntron = 200000; //200kb
$sGenomes = "genomes.txt";
$sSnpEffTablePrefix = "snpeff_ret";
$sTmpDir = "confirmmissing_tmp";
$sOutFile = "confirmmissing_ret.txt";
$sBLASTDIR = "/public/software/conda_envs/funannotate1.8.15/bin/";
$nThreads = 20;
$nAlnPercComplete = 0.9; // if lower than this

$hO = fopen($sOutFile, 'w');

list( $arrOrthList, $arrSampleIDs) = fnLoadOrthList($sSnpEffTablePrefix) ;
print_r($arrSampleIDs);
//print_r($arrOrthList);
//die();
$arrGenomes = fnLoadAllGenomesGFFs($sGenomes);
//die();
fwrite($hO, "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\t".implode("\t",$arrSampleIDs )."\t".implode("\t",$arrSampleIDs )."\n");

$nblastID = 0;

foreach($arrOrthList as $oOrth) {
	$arrSnpEffRet = array();
	$arrMeta = $oOrth['meta'];
	$nOrID = $arrMeta[2];
	list( $sBestSp , $sBestSyn, $sBestProtID) =  explode("|", $oOrth['best']);
	$nBestP = $oOrth['bestP'];
	$sWD = "$sTmpDir/$nOrID";
	exec("mkdir -p $sWD");

	//write down the best gff and and region sequence:
	if (!array_key_exists($sBestProtID,$arrGenomes[$sBestSp]['gff']['rna2gene'] )) {
		die("Error! $sBestProtID not found in $sBestSp GFF\n".$oOrth['best']);
	}

	$sGENEID = $arrGenomes[$sBestSp]['gff']['rna2gene'][$sBestProtID];
	$oGene = $arrGenomes[$sBestSp]['gff']['genes'][$sGENEID];
	$nOrigStart = $oGene['start']; //"scf" => "", "start"=>0, "end"=>0, "orient" => "", "RNAs" => array()
	$nCoordOffset = $nOrigStart - 1;

	$sPrefixScf = "$sBestSp|".$oGene['scf'];
	$sGFF = "$sPrefixScf\tfunannotate\tgene\t".($oGene['start']-$nCoordOffset)."\t".($oGene['end']-$nCoordOffset)."\t.\t".$oGene['orient']."\t.\tID=$sGENEID;\n"; //bahahascfF0_38  funannotate     gene    4219    6936    .       +       .       ID=btpf0000001;

	foreach($oGene['RNAs'] as $oRNA) { //array('def' => $arrF, 'exons' => array(), 'CDS' => array() );
		$oLn = $oRNA['def'];
		$oLn[0] = $sPrefixScf;
		$oLn[3] -= $nCoordOffset ;
		$oLn[4] -= $nCoordOffset ;
		$sGFF .= implode("\t", $oLn)."\n";
		foreach($oRNA['exon'] as $oLn) {
			$oLn[0] = $sPrefixScf;
			$oLn[3] -= $nCoordOffset ;
			$oLn[4] -= $nCoordOffset ;
			$sGFF .= implode("\t", $oLn)."\n";
		}
		foreach($oRNA['CDS'] as $oLn) {
			$oLn[0] = $sPrefixScf;
			$oLn[3] -= $nCoordOffset ;
			$oLn[4] -= $nCoordOffset ;
			$sGFF .= implode("\t", $oLn)."\n";
		}
	}

	$sGFFFile = "$sWD/best.gff";
	file_put_contents($sGFFFile, $sGFF);

	$sFafile = "$sWD/best.fa";
	file_put_contents( $sFafile, ">$sPrefixScf\n" . wordwrap(substr($arrGenomes[$sBestSp]['genome'][$oGene['scf']], $oGene['start']-1, $oGene['end'] - $oGene['start'] + 1 ), 75,  "\n", true ) . "\n" );
	$nBestGeneLenOnScf = $oGene['end'] - $oGene['start'] + 1;

	//build snpeff database
	exec("mkdir -p $sWD/data/genomes; mkdir -p $sWD/data/best; ln -sf `realpath $sFafile`  $sWD/data/genomes/; ln -sf `realpath $sGFFFile` $sWD/data/best/genes.gff");
	$sSNPEffConfig = "$sWD/snpeff.config";
	file_put_contents($sSNPEffConfig, "data.dir = ./data/\n\nbest.genome : best\n" );
	exec("cd $sWD; $sSNPEFF build -c snpeff.config -dataDir `pwd`/data/  -gff3 -noCheckProtein -noCheckCds -v best >/dev/null 2>&1; $sGFF2Protein best.gff best.fa prot > best.prot.fa");

	foreach($oOrth['fields'] as $sSp => $arrInfo) {

		$nAlnPercSyn = fnGetAlnPerc($arrInfo[0]);
		$nAlnPercNS = fnGetAlnPerc($arrInfo[1]);

		if ($nAlnPercSyn >=$nAlnPercComplete || $nAlnPercNS >=$nAlnPercComplete) { //not both missing
			if ($arrInfo[0]!='NA') {
				$arrSnpEffRet[$sSp."|SYN"] = array('raw' => $arrInfo[0]); // copy over
			}

			if ($arrInfo[1]!='NA') {
				$arrSnpEffRet[$sSp."|NS"] = array('raw' => $arrInfo[1]); // copy over
			}
			continue;
		}


		$sIsSyn = "NS";
		$sPrevRunInfo = (($nAlnPercSyn > $nAlnPercNS)? $arrInfo[0]:$arrInfo[1]);
		//try to search for a best alignment
		exec("cd $sWD; $sBLASTDIR/tblastn -num_threads $nThreads -db ../$sSp -query best.prot.fa -out $sSp.txt -outfmt \"7 sseqid ssac qstart qend sstart send pident evalue bitscore\" > screen.log 2>&1; cat $sSp.txt  | grep -v '#' | ". 'awk \'{if ($6>=90) { if ($4<$5) { print $1"\t"$4"\t"$5"\t.\t"$6"\t+";} else {print $1"\t"$5"\t"$4"\t.\t"$6"\t-"} } }\' | sort -k1,1 -k2,2n > '.$sSp.'.candidate.bed; ' . " $sBEDTOOLS merge -i $sSp.candidate.bed -d $nMaxIntron  > $sSp.merged.bed;");

		$oBestHit = fnGetBestBlastHit("$sWD/$sSp.merged.bed");
		//die();

		if ($oBestHit === false) {
			continue;
		}

		$sProtID = $sThisGENEID = "blast".($nblastID++) . "=".$oBestHit[0].":".$oBestHit[1]."-".$oBestHit[2]; //$arrGenomes[$sSp]['gff']['rna2gene'][$sProtID];
		//$oThisGene = $arrGenomes[$sSp]['gff']['genes'][$sThisGENEID];
		$sThisFaFile = "$sWD/other.$sSp.$sIsSyn.$sProtID.fa";

		$nCutStart = $oBestHit[1]-1 - $nBestGeneLenOnScf;
		$nCutEnd = $oBestHit[2]-1 + $nBestGeneLenOnScf;

		$nCutStart = ($nCutStart < 0)? 0:$nCutStart;

		$sIDString = "$sSp|$sIsSyn|$sThisGENEID";

		file_put_contents($sThisFaFile, ">$sIDString\n" .wordwrap(substr($arrGenomes[$sSp]['genome'][$oBestHit[0]], $nCutStart,$nCutEnd- $nCutStart + 1 ), 75,  "\n", true ) . "\n" );
		//die();
		//perform alignment
		exec("cd $sWD; bash ../../aln.sh best.fa 'other.$sSp.$sIsSyn.$sProtID.fa' best '$sSp.$sIsSyn.$sProtID' >/dev/null 2>&1");
		//check if there is alignment:
		if (filesize("$sWD/pairwise.aln.best.$sSp.$sIsSyn.$sProtID.paf") < 10) {
			//no alignment
			echo("$nOrID $sIDString produces no alignment with the best version\n");
			$arrSnpEffRet[$sSp."|".$sIsSyn] = array('aln' => array('alnlen' => 0, 'alnperc' =>0, 'noalignment' => 1) , 'vartypes' => array() , 'effs' => array(), 'protid' => $sProtID, 'lenpval' => $arrPVals[$sSp][$nOrID][$nIsSyn] , 'prevruninfo' => $sPrevRunInfo, 'prevrunalnperc' => max($nAlnPercSyn, $nAlnPercNS));
			continue;
		}

		//check alignment length
		$hPAF = fopen("$sWD/pairwise.aln.best.$sSp.$sIsSyn.$sProtID.paf", 'r');

		$nAlnLen = 0;
		while(false !== ($sLn = fgets($hPAF) ) ) {
			$sLn = trim($sLn);
			if ($sLn == '') continue;
			$arrPAF = explode("\t", $sLn);
			$nAlnLen += $arrPAF[10] ;
			$nTotalLen = $arrPAF[6];
			$nAlnPerc = $nAlnLen/$nTotalLen;
		}

		exec("cd $sWD; $sSNPEFF ann -c snpeff.config best 'out.best.$sSp.$sIsSyn.$sProtID.vcf.gz' | bgzip -c > 'ann.$sSp.$sIsSyn.$sProtID.vcf.gz'");
		//die();
		$hAnn = popen("zcat $sWD/ann.$sSp.$sIsSyn.$sProtID.vcf.gz", 'r');
		$arrVarTypes = array();
		$arrDelEffects = array();
		while(false!== ($sLn = fgets($hAnn) )) {
			$sLn = trim($sLn);
			if ($sLn == '' || $sLn[0] =='#') continue;
			$arrF = explode("\t", $sLn);
			$arrAnn = fnParseAnnot($arrF[7]);
			if (!array_key_exists("ANN", $arrAnn)) {
				continue;
			}
			$arrEffAnn = explode("|", $arrAnn['ANN']);
			$arrT = explode("&", $arrEffAnn[1]);
			$arrE = explode("&", $arrEffAnn[2]);
			foreach($arrT as $nIdx => $sT) {

				if (!array_key_exists($sT,$arrVarTypes  )) {
					$arrVarTypes[$sT] = 0;
				}


				$arrVarTypes[$sT] += 1;

			}

			foreach($arrE as $nIdx => $sE) {
				if (!array_key_exists($sE,$arrDelEffects  )) {
					$arrDelEffects[$sE] = 0;
				}
				$arrDelEffects[$sE] += 1;
			}

		}
		pclose($hAnn);

		$arrSnpEffRet[$sSp."|".$sIsSyn] = array('aln' => array('alnlen' => $nAlnLen, 'alnperc' => $nAlnPerc) , 'vartypes' => $arrVarTypes , 'effs' => $arrDelEffects, 'protid' => $sProtID , 'prevruninfo' => $sPrevRunInfo, 'prevrunalnperc' => max($nAlnPercSyn, $nAlnPercNS)  );
	}

	//print_r($arrSnpEffRet);
	//die();
	fwrite($hO, implode("\t", $arrMeta)."\t" );
	$arrWriteCols = array_fill(0, count($arrSampleIDs)*2, "NA"); 
	foreach($arrSampleIDs as $nIdx => $sSample) {
		foreach(array('SYN', 'NS') as $sSynType ) {
			$sSearch = "$sSample|$sSynType";
			$nSetCol = ($sSynType=='SYN')? $nIdx : $nIdx+count($arrSampleIDs) ;

			if ($sSearch == "$sBestSp". "|". "$sBestSyn") {
				$arrWriteCols[$nSetCol] = "BEST|"."$sBestProtID"."="."$nBestP";
				continue;
			}
			if (!array_key_exists($sSearch , $arrSnpEffRet) ) {
				$arrWriteCols[$nSetCol] = "NA";
				continue;
			}

			$oSnpEff = $arrSnpEffRet[$sSearch];

			if (array_key_exists('raw', $oSnpEff)) {
				$arrWriteCols[$nSetCol] = $oSnpEff['raw']; // copy over previous results
			} else {
				if ($oSnpEff['aln']['alnperc'] > $oSnpEff['prevrunalnperc']) { // if alignment percentage is longer
					$arrOutField = array();
					$arrOutField[] = fnAnn2Str($oSnpEff['aln']); //wrong
					$arrOutField[] = fnAnn2Str($oSnpEff['vartypes']);
					$arrOutField[] = fnAnn2Str($oSnpEff['effs']);
					$arrOutField[] = $oSnpEff['protid'];				
					$arrWriteCols[$nSetCol] = implode("|", $arrOutField);
				} else {
					$arrWriteCols[$nSetCol] = $oSnpEff['prevruninfo']; // copy over previous results
				}
			}
		}
	}

	fwrite($hO, implode("\t", $arrWriteCols)."\n" );


}

function fnGetAlnPerc($sInfo) {

	if (strpos($sInfo, 'BEST|') !== false) {
		return 1;
	} 
	preg_match('/alnperc=([0-9\.]+)/', $sInfo, $arrM);
	if (count($arrM)!=2) {
		return -1;
	}

	return $arrM[1];
}

function fnAnn2Str($arr) {
	$arrRet = array();
	foreach($arr as $key => $val) {
		$arrRet[] = "$key"."=". "$val";
	}

	return implode(";", $arrRet);
}

function fnGetBestBlastHit($s) {
	$h = fopen($s, 'r');
	$arrHits = array();
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$nHitLen = $arrF[2] - $arrF[1];
		$arrHits[$nHitLen] = $arrF;
	}

	if (count($arrHits) == 0) {
		return false;
	}
	return $arrHits[max(array_keys($arrHits))];
}

function fnLoadOrthList($sPrefix) {
	$hPval = fopen("$sPrefix.txt", 'r');

	$nSampleCount = 0;
	$arrSampleIDs = array();
	$arrCol2ID = array();
	
	$arrMeta = array();
	while(false !== ($sLn = fgets($hPval) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == "pgChr") {
			$nSampleCount = (count($arrF) - 4)/2;
			$arrSampleIDs = array_slice($arrF, 4, $nSampleCount);
			foreach($arrSampleIDs as $nIdx => $sID) {
				$arrCol2ID[4+$nIdx ] = $sID;
				//$arrCol2ID[5+$nSampleCount+$nIdx ] = $sID;
			}
			continue;
		}

		if ($arrF[3] == 'NA') {
			continue;
		}
		$sOrID = $arrF[2];
		$oOrth = array();
		$oOrth['meta'] = array_slice($arrF, 0,4);
		$arrPVals = array();
		foreach($arrCol2ID as $nCol => $sID) {
			$nColSupp = $nCol+$nSampleCount;
			$arrPVals[$sID] = array('NA', 'NA');
			if (array_key_exists($nCol, $arrF)) {
				$arrPVals[$sID][0] = $arrF[$nCol];
				$oBest = isBest($arrF[$nCol]);
				if ($oBest!==false) {
					$oOrth['best'] = "$sID|SYN|".$oBest[0];
					$oOrth['bestP'] = $oBest[1];
				}
			}
			if (array_key_exists($nColSupp, $arrF)) {
				$arrPVals[$sID][1] = $arrF[$nColSupp];
				$oBest = isBest($arrF[$nColSupp]);
				if ($oBest!==false) {
					$oOrth['best'] = "$sID|NS|".$oBest[0];
					$oOrth['bestP'] = $oBest[1];
					
				}
			}

		}
		$oOrth['fields'] = $arrPVals;
		$arrMeta[$sOrID] = $oOrth;
	}


	return array( $arrMeta , $arrSampleIDs);


}

function isBest($s) {
	$arrF = explode('|', $s);
	if ($arrF[0] == 'BEST') {
		if (count($arrF)<2) {
			die($s);
		}
		return explode("=", $arrF[1]);
	}

	return false;
}

function fnLoadAllGenomesGFFs($s) {
	global $sTmpDir, $sBLASTDIR;
	$arrRet = array();

	if (!file_exists($s)) {
		die("file $s not found\n");
	}
	$h = fopen($s, 'r');

	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn );

		if ($sLn == '') continue;

		$arrF = explode("\t", $sLn);
		if (count($arrF)!=3) {
			die("Format error: \n$sLn\n");
		}
		$arrRet[$arrF[0]] = array('genome' => readFasta($arrF[1]), 'gff' => readGFFGenes($arrF[2]) ); 
		if (!file_exists($sTmpDir."/".$arrF[0].".ndb")) {
			exec("FA=`realpath ".$arrF[1]."`; mkdir -p $sTmpDir; cd $sTmpDir; $sBLASTDIR/makeblastdb -in \$FA -dbtype nucl -out ".$arrF[0].";");
		}
	}

	return $arrRet;
}
function readGFFGenes($s) {
	echo("reading $s ...\n");
	if (!file_exists($s)) {
		die("file $s not found\n");
	}

	$h = fopen($s, 'r');
	$arrRet = array();
	$arrRNA2Gene = array();
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn );

		if ($sLn == '') continue;

		$arrF = explode("\t", $sLn);
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseAnnot($arrF[8]);
			if (!array_key_exists("ID", $arrAnnot) || !array_key_exists("Parent", $arrAnnot) ) {
				continue;
			}
			$sRNAID = $arrAnnot['ID'];
			$sGENEID = $arrAnnot['Parent'];
			if (!array_key_exists($sGENEID, $arrRet)) {
				$arrRet[$sGENEID] = array("scf" => "", "start"=>0, "end"=>0, "orient" => "", "RNAs" => array() );
			}

			if ($arrRet[$sGENEID]['scf'] == "") {
				$arrRet[$sGENEID]['scf'] = $arrF[0];
				$arrRet[$sGENEID]['orient'] = $arrF[6];
				$arrRet[$sGENEID]['start'] = $arrF[3];
				$arrRet[$sGENEID]['end'] = $arrF[4];
			}

			$arrRet[$sGENEID]['start'] = min( $arrF[3] , $arrRet[$sGENEID]['start']) ;
			$arrRet[$sGENEID]['end'] = max($arrF[4] , $arrRet[$sGENEID]['end']) ;

			$arrRet[$sGENEID]['RNAs'][$sRNAID] = array('def' => $arrF, 'exons' => array(), 'CDS' => array() );
			$arrRNA2Gene[$sRNAID] = $sGENEID;
		}

		if ($arrF[2] == 'CDS' || $arrF[2] == 'exon') {
			$arrAnnot = fnParseAnnot($arrF[8]);
			$sRNAID = $arrAnnot['Parent'];
			if (array_key_exists($sRNAID, $arrRNA2Gene)) {
				$sGENEID = $arrRNA2Gene[$sRNAID];
				$arrRet[$sGENEID]['RNAs'][$sRNAID][$arrF[2]][] = $arrF;
			}
		}

	}

	return array( 'genes' => $arrRet , 'rna2gene' => $arrRNA2Gene);
}

function readFasta($s) {
		echo("reading $s ...\n");
		if (!file_exists($s)) {
			die("file $s not found\n");
		}
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
                        $sSeq = '';
                        continue;
                }

                $sSeq .= $sLn;

        } while(true);

        return $arrRet;
}

function fnParseAnnot($s) {
	$arrF = explode(";", $s);
	$arrRet = array();
	foreach($arrF as $sF) {
		if (trim($sF) == '') continue;
		list($sName, $sVal) = explode("=",$sF, 2);
		$arrRet[trim($sName, " '\"")] = trim($sVal, " \t'\"");
	}

	return $arrRet;
}

?>
