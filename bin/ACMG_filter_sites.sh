inputvcf=$1
Manuel_evidence=$2
outdir=$3


Script_dir=$(dirname "$(readlink -f "$0")")
input_dir=$(dirname "$(readlink -f "$Manuel_evidence")")

cd $outdir/result/
python3.6 $Script_dir/ACMG_filter.py $inputvcf ACMG_sites_af.tsv
awk -F'\t' '{print $0}' ACMG_sites_af.tsv | cut -f1-4 >sitelist
cat $outdir/samplelabling/*.classfy.tsv | head -n 1 | sed 's/^/sampleid\t/' > classfy_results_samples.temp0
for i in `ls $outdir/samplelabling/*.tsv`
do
    sample=`basename $i|sed 's/.classfy.tsv//'`
    grep -v "chrom" $i | sed 's/^/'$sample'\t/g' >> classfy_results_samples.temp0
done
cut -f1 $input_dir/transid.txt > transcript_list.txt
grep -wFf transcript_list.txt classfy_results_samples.temp0 >classfy_results_samples.temp1
grep -wFf sitelist classfy_results_samples.temp1 > classfy_results_samples.temp2
awk -F'\t' 'NR==FNR{a[$1$2$3$4$5]=$7"\t"$8"\t"$9"\t"$10"\t"$11}NR>FNR{if(a[$2$3$4$5$6]){print $0"\t"a[$2$3$4$5$6]}else{print $0"\tNA\tNA\tNA\tNA\tNA"}}' ACMG_split_feature.tsv classfy_results_samples.temp2 > classfy_results_samples.temp3
cat $outdir/samplelabling/*.classfy.tsv | head -n 1 | sed 's/$/\thgvsc\thgvsp\timpact\tconsequence\tvep_splice_region/'| sed 's/^/sampleid\t/' | cat - classfy_results_samples.temp3 > ACMG_classfy_results_filtered.tsv
head -1 ACMG_classfy_results_filtered.tsv > ACMG_classfy_results_filtered.1.tsv
awk -F'\t' '$17!="./0"&&$17!="0/0"&&$17!="0/." {print $0}' ACMG_classfy_results_filtered.tsv|sed '1d' |sort -u >> ACMG_classfy_results_filtered.1.tsv


cat ACMG_classfy_results_filtered.1.tsv | awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR&&FNR==1{for(i=6;i<=NF;i++){ix[i]=$i}}NR==FNR&&NR>=2{for(i=6;i<=NF;i++){if($i!="NA"){a[$1$2$3$4$5]=a[$1$2$3$4$5]""";"""ix[i]"="$i}}}NR>FNR{gsub(/PP5=[0-9]+;?/, "", $29);gsub(/PP5=[0-9]+;?/, "", $25);if($(NF-1)!~/missense/&&$25~/PS1/&&$(NF-1)!~/splice/){gsub(/PS1=[0-9]+;?/, "", $25);gsub(/PS1=[0-9]+;?/, "", $29)};if($(NF-1)~/synonymous/&&$25!~/PP3/){$25=$25";BP7=3";$29=$29";BP7=3"};if(substr($25,length($25))==";"){$25=substr($25,1,length($25)-1)};if(substr($29,length($29))==";"){$29=substr($29,1,length($29)-1)};if(substr($25,1,1)==";"){$25=substr($25,2)};if(substr($29,1,1)==";"){$29=substr($29,2)};gsub(/;;/, ";", $25);gsub(/;;/, ";", $29);if($1$2$3$4$5 in a){$29=$29""a[$1$2$3$4$5];$25=$25""a[$1$2$3$4$5];print $0}else{print $0}}' $Manuel_evidence - >ACMG_classfy_results_filtered.final.tsv


#cat ACMG_classfy_results_filtered.1.tsv | awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR&&FNR==1{for(i=6;i<=NF;i++){ix[i]=$i}}NR==FNR&&NR>=2{for(i=6;i<=NF;i++){if($i!="NA"){a[$1$2$3$4$5]=a[$1$2$3$4$5]""";"""ix[i]"="$i}}}NR>FNR{gsub(/PP7=[0-9]+;?/, "", $29);gsub(/PP5=[0-9]+;?/, "", $29);gsub(/PP7=[0-9]+;?/, "", $25);gsub(/PP5=[0-9]+;?/, "", $25);if($(NF-1)!~/missense/&&$25~/PS1=/&&$(NF-1)!~/splice/){gsub(/PS1=[0-9]+;?/, "", $25);gsub(/PS1=[0-9]+;?/, "", $29)};if($(NF-1)~/synonymous/&&$25!~/PP3/){$25=$25";BP7=3";$29=$29";BP7=3"};if(substr($25,length($25))==";"){$25=substr($25,1,length($25)-1)};if(substr($29,length($29))==";"){$29=substr($29,1,length($29)-1)};if($1$2$3$4$5 in a){$29=$29""a[$1$2$3$4$5];$25=$25""a[$1$2$3$4$5];print $0}else{print $0}}' $Manuel_evidence - >ACMG_classfy_results_filtered.final.tsv


Rscript $Script_dir/Classification.R ACMG_classfy_results_filtered.final.tsv Classification.tsv
cat Classification.tsv | cut -f1-7,10,11,17,19,20,25,29,30-40 >FIA_classification.tsv

rm classfy_results_*
rm Classification.tsv
rm ACMG_classfy_results_filtered.*
rm ../none_sample*
rm ../all_prediction_20181105.txt
rm -rf ../data
rm -rf ../nonsampledir_both
rm -rf ../nonsample_evidence_ad
rm -rf ../nonsample_evidence_ar
rm -rf ../nonsample_evidence_other
rm -rf ../P_satat
rm -rf ../sample_evidence
