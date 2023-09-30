
FALCONPHASEBIN=/data/software/FALCON-Phase/bin
Pcontigs=../falcon/4-polish/cns-output/polished_p_ctgs.fasta
Hcontigs=../falcon/4-polish/cns-output/polished_h_ctgs.fasta
#Pcontigs=../falcon/3-unzip/all_p_ctg.fasta
#Hcontigs=../falcon/3-unzip/all_h_ctg.fasta

#python3 $FALCONPHASEBIN/preprocess_diploid_asm_for_fc_phase.py -P $Pcontigs -H $Hcontigs -o cleaned_contigs > omitted_htigs.txt 2>&1 # writes name_mappings.txt
#please edit config.json before running !
bash run.sh > falconphasepipeline.log 2>&1
