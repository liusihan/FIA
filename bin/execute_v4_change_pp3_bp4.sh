#!/usr/bin/bash


input1=$1

input2=$2


ACMGDir_AD=$3
ACMGDir_AR=$4


script=$0
scriptdir=`dirname $script`
ACMGDir=$ACMGDir_AR/temp
mkdir -p $ACMGDir

cd $ACMGDir


output1=SNV_allsites_v1.rmsk.tsv

bcftools query \
-H \
-f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/gnomAD_AF_SAS_WES\t%INFO/gnomAD_grpmax_WES\t%INFO/gnomAD_grpmax_WGS\t%INFO/gnomAD_AF_EAS_WES\t%INFO/gnomAD_AF_ASJ_WES\t%INFO/gnomAD_AF_MID_WES\t%INFO/gnomAD_AF_FIN_WES\t%INFO/gnomAD_AF_NFE_WES\t%INFO/gnomAD_AF_AFR_WES\t%INFO/DBNSFP_INTERPRO_DOMAIN\t%INFO/VEP_MaxEntScan_diff\t%INFO/VEP_Consequence\t%INFO/VEP_DOMAINS\t%INFO/VEP_Feature\t%INFO/VEP_IMPACT\t%INFO/VEP_HGVSc\t%INFO/VEP_HGVSp\t%INFO/VEP_PICK\t%INFO/VEP_Gene\t%INFO/gnomAD_AF_AMI_WES\t%INFO/clinvar_CLNSIG\t%INFO/DVD_PATHOGENICITY\t%INFO/DBNSFP_REVEL_SCORE\t%INFO/rmsk_repName\n' \
$input1 \
> \
$output1

echo "Done: query vcf-->tsv"


perl $scriptdir/transcript_split.pl $ACMGDir

perl $scriptdir/ACMG_tag_AR_20190812.pl $ACMGDir $input2

perl $scriptdir/ACMG_tag_AD_20190812.pl $ACMGDir $input2

perl $scriptdir/sed_v2.pl $ACMGDir_AD $ACMGDir_AR
