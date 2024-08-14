inputvcf=$1
outdir=$2

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DBNSFP_REVEL_SCORE\t%INFO/VEP_MaxEntScan_diff\n' -H $inputvcf > $outdir/data/selected.tsv

script=$0

scriptdir=`dirname $script`

cat $outdir/data/selected.tsv | awk 'BEGIN{FS="\t";OFS="\t";print "chrom","pos","ref","alt","PP3","BP4"}NR>=FNR{split($6,a,"|");$6=a[2];if($5>=0.7||$6>=6){print $1,$2,$3,$4,3,"NA"}else if(($5!="."&&$5<=0.15)||($6!=""&&$6<6)){print $1,$2,$3,$4,"NA",3}else{print $1,$2,$3,$4,"NA","NA"}}' > $outdir/data/PP3_BP4.tsv

python3.6 $scriptdir/ACMG_split_feature.py $inputvcf $outdir/result/ACMG_split_feature.tsv

head -1 $outdir/data/PP3_BP4.tsv | sed 's/alt/alt\ttranscript/' >$outdir/nonsample_evidence_other/PP3_BP4.tsv
awk -F'\t' 'NR==FNR{a[$1$2$3$4]=$5"\t"$6}NR>FNR{if(a[$1$2$3$4]){OFS="\t";print $1,$2,$3,$4,$5,a[$1$2$3$4]}}' $outdir/data/PP3_BP4.tsv $outdir/result/ACMG_split_feature.tsv >> $outdir/nonsample_evidence_other/PP3_BP4.tsv

