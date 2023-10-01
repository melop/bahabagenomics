<?php
ini_set('memory_limit', '50000M');

$sSNPEFF = "/public/software/snpEff/exec/snpeff";
$sGenomes = "genomes.txt";
$sOrthOutlierTablePrefix = "orth_length_outliers";
$sTmpDir = "snpeff_tmp";
$sOutFile = "snpeff_ret.txt";

$hO = fopen($sOutFile, 'w');

list( $arrOrthList, $arrSampleIDs, $arrPVals ) = fnLoadOrthList($sOrthOutlierTablePrefix) ;
$arrGenomes = fnLoadAllGenomesGFFs($sGenomes);
fwrite($hO, "pgChr\tpgOrd\tpgID\tOutgroupAlnPercMean\t".implode("\t",$arrSampleIDs )."\t".implode("\t",$arrSampleIDs )."\n");

foreach($arrOrthList as $oOrth) {
	$arrSnpEffRet = array();
	$arrMeta = $oOrth['meta'];
	echo("Processing: ". implode(" ", $arrMeta)."\n");
	$nOrID = $arrMeta[2];
	list( $sBestSp , $sBestSyn, $sBestProtID) =  explode("|", $oOrth['best']);
	$nBestP = $oOrth['bestP'];
	$sWD = "$sTmpDir/$nOrID";
	exec("mkdir -p $sWD");

	//write down the best gff and and region sequence:
	if (!array_key_exists($sBestProtID,$arrGenomes[$sBestSp]['gff']['rna2gene'] )) {
		die("Error! $sBestProtID not found in $sBestSp GFF\n");
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
	exec("cd $sWD; $sSNPEFF build -c snpeff.config -dataDir `pwd`/data/  -gff3 -noCheckProtein -noCheckCds -v best >/dev/null 2>&1");

	foreach($oOrth['other'] as $sIDString => $nP) {

		list($sSp, $sIsSyn, $sProtID) = explode("|", $sIDString) ;

		$nIsSyn = ($sIsSyn=='SYN')? 0:1;

		/* look at identical ones anyways
		if ($nP == $nBestP) {
			echo("$nOrID $sIDString Protein appears identical with the best version\n");
			$arrSnpEffRet[$sSp."|".$sIsSyn] = array('aln' => array('alnlen' => 1, 'alnperc' =>1, 'identical' => 1) ,  'vartypes' => array() , 'effs' => array(), 'protid' => $sProtID, 'lenpval' => $arrPVals[$sSp][$nOrID][$nIsSyn] );
			continue;
		}
		*/

		$sThisGENEID = $arrGenomes[$sSp]['gff']['rna2gene'][$sProtID];
		$oThisGene = $arrGenomes[$sSp]['gff']['genes'][$sThisGENEID];
		$sThisFaFile = "$sWD/other.$sSp.$sIsSyn.$sProtID.fa";

		$nCutStart = $oThisGene['start']-1 - $nBestGeneLenOnScf;
		$nCutEnd = $oThisGene['end']-1 + $nBestGeneLenOnScf;

		$nCutStart = ($nCutStart < 0)? 0:$nCutStart;

		file_put_contents($sThisFaFile, ">$sIDString\n" .wordwrap(substr($arrGenomes[$sSp]['genome'][$oThisGene['scf']], $nCutStart,$nCutEnd- $nCutStart + 1 ), 75,  "\n", true ) . "\n" );
		//die();
		//perform alignment
		exec("cd $sWD; bash ../../aln.sh best.fa 'other.$sSp.$sIsSyn.$sProtID.fa' best '$sSp.$sIsSyn.$sProtID' > aln.$sSp.$sIsSyn.$sProtID.log 2>&1");
		//check if there is alignment:
		if (filesize("$sWD/pairwise.aln.best.$sSp.$sIsSyn.$sProtID.paf") <= 0) {
			//no alignment
			echo("$nOrID $sIDString produces no alignment with the best version\n");
			$arrSnpEffRet[$sSp."|".$sIsSyn] = array('aln' => array('alnlen' => 0, 'alnperc' =>0, 'noalignment' => 1) , 'vartypes' => array() , 'effs' => array(), 'protid' => $sProtID, 'lenpval' => $arrPVals[$sSp][$nOrID][$nIsSyn] );
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

		$arrSnpEffRet[$sSp."|".$sIsSyn] = array('aln' => array('alnlen' => $nAlnLen, 'alnperc' => $nAlnPerc) , 'vartypes' => $arrVarTypes , 'effs' => $arrDelEffects, 'protid' => $sProtID, 'lenpval' => $arrPVals[$sSp][$nOrID][$nIsSyn] );
	}

	//print_r($arrSnpEffRet);
	//die();
	fwrite($hO, implode("\t", $arrMeta)."\t" );
	$arrWriteCols = array_fill(0, count($arrSampleIDs)*2, "NA");
	foreach($arrSampleIDs as $nIdx => $sSample) {
		foreach(array('SYN', 'NS') as $sSynType ) {
			$sSearch = "$sSample|$sSynType";
			$nSetCol = ($sSynType=='SYN')? $nIdx : $nIdx+count($arrSampleIDs);

			if ($sSearch == "$sBestSp". "|". "$sBestSyn") {
				$arrWriteCols[$nSetCol] = "BEST|"."$sBestProtID"."="."$nBestP";
				continue;
			}
			if (!array_key_exists($sSearch , $arrSnpEffRet) ) {
				$arrWriteCols[$nSetCol] = "NA";
				continue;
			}

			$oSnpEff = $arrSnpEffRet[$sSearch];

			$arrOutField = array();
			$arrOutField[] = fnAnn2Str($oSnpEff['aln']); //wrong
			$arrOutField[] = fnAnn2Str($oSnpEff['vartypes']);
			$arrOutField[] = fnAnn2Str($oSnpEff['effs']);
			$arrOutField[] = $oSnpEff['protid']."=".$oSnpEff['lenpval'];

			$arrWriteCols[$nSetCol] = implode("|", $arrOutField);
		}
	}

	fwrite($hO, implode("\t", $arrWriteCols)."\n" );


}

function fnAnn2Str($arr) {
	$arrRet = array();
	foreach($arr as $key => $val) {
		$arrRet[] = "$key"."=". "$val";
	}

	return implode(";", $arrRet);
}

function fnLoadOrthList($sPrefix) {
	$hPval = fopen("$sPrefix.txt", 'r');
	$hIDs = fopen("$sPrefix.ids.txt", 'r');

	$nSampleCount = 0;
	$arrSampleIDs = array();
	$arrCol2ID = array();
	$arrPVals = array();
	while(false !== ($sLn = fgets($hPval) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == "pgChr") {
			$nSampleCount = (count($arrF) - 5)/2;
			$arrSampleIDs = array_slice($arrF, 5, $nSampleCount);
			foreach($arrSampleIDs as $nIdx => $sID) {
				$arrPVals[$sID] = array();
				$arrCol2ID[5+$nIdx ] = $sID;
				//$arrCol2ID[5+$nSampleCount+$nIdx ] = $sID;
			}
			continue;
		}

		if ($arrF[3] == 'NA') {
			continue;
		}
		$sOrID = $arrF[2];
		foreach($arrCol2ID as $nCol => $sID) {
			$nColSupp = $nCol+$nSampleCount;
			$arrPVals[$sID][$sOrID] = array('NA', 'NA');
			if (array_key_exists($nCol, $arrF)) {
				$arrPVals[$sID][$sOrID][0] = $arrF[$nCol];
			}
			if (array_key_exists($nColSupp, $arrF)) {
				$arrPVals[$sID][$sOrID][1] = $arrF[$nColSupp];
			}

		}
	}

	$arrRet = array();
	while(false !== ($sLn = fgets($hIDs) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == "pgChr") {
			foreach($arrCol2ID as $nCol => $sID) {
				if ($arrF[$nCol] != $sID || $arrF[$nCol+$nSampleCount]!=$sID) {
					die("Error: ID mismatches in header of $sPrefix\n");
				}
			}
			continue;
		}

		$sOrID = $arrF[2];

		$oRet = array('meta' => array_slice($arrF,0,4), 'best' => '', 'bestP' => 0, 'other'=>array());

		$arrAllGenes = array();
		foreach($arrCol2ID as $nCol => $sID) {
			$sSynID = $arrF[$nCol];
			$sNSID = $arrF[$nCol+$nSampleCount];
			$nSynP = $arrPVals[$sID][$sOrID][0];
			$nNSP = $arrPVals[$sID][$sOrID][1];
			if ($sSynID != 'NA' && $nSynP != 'NA') {
				$arrAllGenes["$sID"."|SYN|".$sSynID] = $nSynP;
			}
			if ($sNSID != 'NA' && $nNSP != 'NA') {
				$arrAllGenes["$sID"."|NS|".$sNSID] = $nNSP;
			}

		}

		if (count($arrAllGenes) == 0) {
			echo("Warning: all genes missing:\t$sLn\n");
			continue;		
		}

		$nMaxP = max($arrAllGenes);

		if ($nMaxP<0.01) {
			echo("Warning: all genes too short or missing:\t$sLn\n");
			continue;
		}

		//max p is set as the reference
		$sBestModel = "";
		foreach($arrAllGenes as $sIDString => $nP) {
			if ($nP == $nMaxP) {
				$sBestModel = $sIDString;
				break;
			}
		}

		$oRet['best'] = $sBestModel;
		$oRet['bestP'] = $nMaxP;

		foreach($arrAllGenes as $sIDString => $nP) {
			if ( $sIDString == $sBestModel) {
				continue;
			}

			$oRet['other'][$sIDString] = $nP;
		}


		$arrRet[] = $oRet;
	}

	return array( $arrRet , $arrSampleIDs,$arrPVals);


}

function fnLoadAllGenomesGFFs($s) {
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
