
paste /home/hjb/FIA_github/test/FIA/nonsampledir_both/*.tsv | cut -f1-6,12-13,19 > /home/hjb/FIA_github/test/FIA/none_sample_evidence_both.tsv
paste /home/hjb/FIA_github/test/FIA/nonsample_evidence_ad/*.tsv > /home/hjb/FIA_github/test/FIA/none_sample_evidence_AD.tsv
paste /home/hjb/FIA_github/test/FIA/nonsample_evidence_ar/*.tsv > /home/hjb/FIA_github/test/FIA/none_sample_evidence_AR.tsv
echo "step2 finished" >> /home/hjb/FIA_github/test/FIA/bin/record.txt

