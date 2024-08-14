
## 

input=$1
outdir=$2
psites=$3
linkagesites=$4
hl_samplelist=$5
ARgenes=$6

script=$0
scriptdir=`dirname $script`




paral_num=20

mkdir -p $outdir/temp



bcftools view -h $input |tail --lines=1|awk '{for(i=10;i<=NF;i++)printf("%s\n",$i);print}'>$outdir/temp/all_samples.txt
cat $outdir/temp/all_samples.txt $hl_samplelist |sort|uniq -d > $outdir/temp/HL_samples.txt
bcftools view --threads $paral_num --force-samples -S $outdir/temp/HL_samples.txt $input -Oz -o $outdir/temp/HL_samples.vcf.gz
tabix -f -p vcf $outdir/temp/HL_samples.vcf.gz

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_Feature\t%INFO/VEP_IMPACT\t%INFO/VEP_SYMBOL\t%INFO/gnomAD_AF_AMI_WES\t%INFO/gnomAD_AF_SAS_WES\t%INFO/gnomAD_grpmax_WES\t%INFO/gnomAD_grpmax_WGS\t%INFO/gnomAD_AF_EAS_WES\t%INFO/gnomAD_AF_ASJ_WES\t%INFO/gnomAD_AF_MID_WES\t%INFO/gnomAD_AF_FIN_WES\t%INFO/gnomAD_AF_NFE_WES\t%INFO/gnomAD_AF_SAS_WES\t%INFO/VEP_SpliceRegion\t%INFO/VEP_Consequence[\t%GT]\n' $outdir/temp/HL_samples.vcf.gz > $outdir/temp/HL_samples_extarct.tsv

python3.6 $scriptdir/pm3_extract_info.py $outdir/temp/HL_samples_extarct.tsv $psites $linkagesites $outdir/temp/network.tsv $ARgenes
cut -f5 $outdir/temp/network.tsv | sort -u > $outdir/temp/uniq_transcripts.txt

awk 'BEGIN{OFS="\t";print "chrom","pos","ref","alt","transcript","PM3"}' > $outdir/pm3.tsv
cat $outdir/temp/uniq_transcripts.txt | xargs -P 20 -n 1 python3.6 $scriptdir/pm3_combine_mul.py $outdir/temp/network.tsv $psites $outdir/pm3.tsv


















