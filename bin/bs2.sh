
## 

input=$1
outdir=$2
psites=$3
samplelist=$4

paral_num=20

mkdir -p $outdir/BS2temp
script=$0
scriptdir=`dirname $script`

## 1. extract nonHL sample
cut -f1 $samplelist >$outdir/BS2temp/control_id.txt
bcftools view --threads $paral_num -S $outdir/BS2temp/control_id.txt --force-samples $input -Oz > $outdir/BS2temp/nonHL_samples.vcf.gz
tabix -f -p vcf $outdir/BS2temp/nonHL_samples.vcf.gz

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_Feature\t%INFO/VEP_IMPACT\t%INFO/VEP_SYMBOL\t%INFO/VEP_Gene\t%INFO/gnomAD_AF_SAS_WES\t%INFO/gnomAD_grpmax_WES\t%INFO/gnomAD_grpmax_WGS\t%INFO/gnomAD_AF_EAS_WES\t%INFO/gnomAD_AF_ASJ_WES\t%INFO/gnomAD_AF_MID_WES\t%INFO/gnomAD_AF_FIN_WES\t%INFO/gnomAD_AF_NFE_WES\t%INFO/gnomAD_AF_AMI_WES[\t%GT]\n' $outdir/BS2temp/nonHL_samples.vcf.gz > $outdir/BS2temp/nonHL_samples_extarct.tsv
## 2. split transcript
python3.6 $scriptdir/extract_info_bs2.py $outdir/BS2temp/nonHL_samples_extarct.tsv $outdir/BS2temp/network.tsv

cut -f5 $outdir/BS2temp/network.tsv | sort -u > $outdir/BS2temp/uniq_transcripts.txt
## 3. network
echo -e "chrom\tpos\tref\talt\ttranscript\tBS2" >$outdir/BS2.tsv
cat $outdir/BS2temp/uniq_transcripts.txt | xargs -P 20 -n 1 python3.6 $scriptdir/bs2_combine_mul_20190717new.py $outdir/BS2temp/network.tsv $psites $outdir/BS2.tsv


















