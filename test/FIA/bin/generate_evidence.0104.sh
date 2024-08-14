
mkdir -p /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/
mkdir -p /home/hjb/FIA_github/test/FIA/sample_evidence/pp4Result/
cp /home/hjb/FIA_github/inputfile/geneset-*.txt /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/
cp /home/hjb/FIA_github/inputfile/phenotype.txt /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/

sh /home/hjb/FIA_github/bin/pp4/getHL_122genes.sh /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/ /home/hjb/FIA_github/test/annotated_raw.AF_stats.vcf.gz
perl /home/hjb/FIA_github/bin/pp4/reformat.pl /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/ HLsamples_122genes.tsv HLsamples_122genes.reformat.tsv
python3.6 /home/hjb/FIA_github/bin/pp4/pp4_main.py /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/phenotype.txt /home/hjb/FIA_github/test/FIA/sample_evidence/rawData/geneset-*.txt "/home/hjb/FIA_github/test/FIA/sample_evidence/pp4Result/" "/home/hjb/FIA_github/test/FIA/sample_evidence/rawData/HLsamples_122genes.reformat.tsv" "/home/hjb/FIA_github/test/annotated_raw.AF_stats.vcf.gz"
echo "step1.4 finished" >> /home/hjb/FIA_github/test/FIA/bin/record.txt

