<?php
//convert from gff3 file to a table listing isoform IDs for each gene ID
$sGFF = "../../../annotations/funannotate1.8.15/hifiasm_hap_F0/btpf0.full.gff3";
$sOut = "BahahaTaipingensisF0.isoforms.txt";

$sGFF = "../../../annotations/funannotate1.8.15/hifiasm_hap_F1/btpf1.full.gff3";
$sOut = "BahahaTaipingensisF1.isoforms.txt";


$sGFF = "../../../annotations/funannotate1.8.15/hifiasm_hap_M0/btpm0.full.gff3";
$sOut = "BahahaTaipingensisM0.isoforms.txt";

$sGFF = "../../../annotations/funannotate1.8.15/hifiasm_hap_M1/btpm1.full.gff3";
$sOut = "BahahaTaipingensisM1.isoforms.txt";


$h = fopen($sGFF, 'r');
$hO = fopen($sOut, 'w');

$arrGenes = array();
while(false!==($sLn = fgets($h) )) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (count($arrF) < 9 ) {
		continue;
	}
	if ($arrF[2] != 'mRNA') {
		continue;
	}

	$arrAnn = fnParseAnnotation($arrF[8]);
	$sID = $arrAnn['ID'];
	$sParent = $arrAnn['Parent'];
	if (!array_key_exists($sParent, $arrGenes)) {
		$arrGenes[$sParent] = array();
	}

	$arrGenes[$sParent][$sID] = true;
}

foreach($arrGenes as $sGeneID => $arrTrans) {
	fwrite($hO, "$sGeneID\t".implode("\t", array_keys($arrTrans))."\n");
}

function fnParseAnnotation($s) {
                $arrMap = array();
                $arrF = explode(';' , $s);
                foreach($arrF as $v) {
                        $arrPair = explode('=' , $v);
                        if (count($arrPair) !=2) continue;
                        $sKey = trim($arrPair[0]);
                        $sValue = preg_replace('/^"|"$/', '', trim($arrPair[1]));
                        $arrMap[$sKey] = $sValue;
                }
                return $arrMap;
        }


?>
