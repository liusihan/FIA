
python /home/hjb/FIA_github/Scripts/FIA_PVS1.py -i /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz -g GRCh37 -o /home/hjb/FIA_github/test/FIA/nonsampledir_both/
perl /home/hjb/FIA_github/bin/addPS1_PM5_by_transcript_20181011.pl -anno /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz -t /home/hjb/FIA_github/inputfile/PUB_P.txt -o /home/hjb/FIA_github/test/FIA/nonsampledir_both/PS1_PM5.tsv
perl /home/hjb/FIA_github/bin/addBP2_by_transcript_20181011.pl -anno /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz -t /home/hjb/FIA_github/inputfile/PUB_P.txt -l /home/hjb/FIA_github/inputfile/BP2_linkage.txt -o /home/hjb/FIA_github/test/FIA/nonsampledir_both/BP2.tsv
##PP3,BP4
sh /home/hjb/FIA_github/bin/pp3_bp4.sh /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz /home/hjb/FIA_github/test/FIA

##
sh /home/hjb/FIA_github/bin/execute_v4_change_pp3_bp4.sh /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz /home/hjb/FIA_github/inputfile/PUB_P.txt /home/hjb/FIA_github/test/FIA/nonsample_evidence_ad /home/hjb/FIA_github/test/FIA/nonsample_evidence_ar
##PM1
sh /home/hjb/FIA_github/bin/pm1.sh /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz /home/hjb/FIA_github/test/FIA/nonsample_evidence_other /home/hjb/FIA_github/inputfile/PM1_hot_region.bed

echo "step1.1 finished" >> /home/hjb/FIA_github/test/FIA/bin/record.txt

