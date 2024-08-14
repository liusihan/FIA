
paste /home/hjb/shliu/ACMG/nonsampledir_both/*.tsv | cut -f1-6,12-13,19 > /home/hjb/shliu/ACMG/none_sample_evidence_both.tsv
paste /home/hjb/shliu/ACMG/nonsample_evidence_ad/*.tsv > /home/hjb/shliu/ACMG/none_sample_evidence_AD.tsv
paste /home/hjb/shliu/ACMG/nonsample_evidence_ar/*.tsv > /home/hjb/shliu/ACMG/none_sample_evidence_AR.tsv
echo "step2 finished" >> /home/hjb/shliu/ACMG/bin/record.txt

