cp /home/hjb/FIA_github/inputfile/all_prediction_20181105.txt /home/hjb/FIA_github/test/FIA/
R --vanilla --no-echo -f /home/hjb/FIA_github/bin/combine_labels_20181207.R --args /home/hjb/FIA_github/test/FIA/none_sample_evidence_both.tsv /home/hjb/FIA_github/test/FIA/none_sample_evidence_AR.tsv /home/hjb/FIA_github/test/FIA/none_sample_evidence_AD.tsv /home/hjb/FIA_github/test/FIA/nonsample_evidence_other
ls /home/hjb/FIA_github/test/FIA/sample_evidence/pp4Result/*.tsv > /home/hjb/FIA_github/test/FIA/samplelist.txt
#mv /home/hjb/FIA_github/test/FIA/ /NAS/shliu/
#ln -s /NAS/shliu/test/FIA /home/hjb/FIA_github/test/FIA 
perl /home/hjb/FIA_github/bin/intervar_samplelist_20190625_statspecific.pl /home/hjb/FIA_github/test/FIA/data/test_dat_ar.txt /home/hjb/FIA_github/test/FIA/data/test_dat_ad.txt /home/hjb/FIA_github/test/FIA/samplelist.txt /home/hjb/FIA_github/inputfile/Manual_evidence.txt /home/hjb/FIA_github/test/FIA/samplelabling/
## variants filtering 
sh /home/hjb/FIA_github/bin/ACMG_filter_sites.sh /home/hjb/FIA_github/test/FIA/data/nonsample.vcf.gz /home/hjb/FIA_github/inputfile/Manual_evidence.txt /home/hjb/FIA_github/test/FIA
### stat
#Rscript /home/hjb/FIA_github/bin/evidence_split.R /home/hjb/FIA_github/test/FIA/result/ /home/hjb/FIA_github/test/FIA/result/ACMG_classfy_results_filtered.forstat.tsv /home/hjb/FIA_github/test/FIA/result/ACMG_split_feature.tsv
#mkdir -p /home/hjb/FIA_github/test/FIA/result/stat/
#Rscript /home/hjb/FIA_github/bin/stat.R /home/hjb/FIA_github/test/FIA/result/stat/ /home/hjb/FIA_github/test/FIA/result/ /public/home/swgenetics_4/database/InputFiles/P_LP_SNP.txt

