#!/usr/bin/bash

outputdir=$1

input=$2

output=$outputdir/HLsamples_122genes.tsv

paral_num=20

cd $outputdir

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_Feature\t%INFO/VEP_SYMBOL\t%INFO/VEP_Gene\t%INFO/gnomAD_AF_AFR_WES\t%INFO/gnomAD_grpmax_WES\t%INFO/gnomAD_grpmax_WGS\t%INFO/gnomAD_AF_EAS_WES\t%INFO/gnomAD_AF_ASJ_WES\t%INFO/gnomAD_AF_MID_WES\t%INFO/gnomAD_AF_FIN_WES\t%INFO/gnomAD_AF_NFE_WES[\t%GT;%AD;%DP;%GQ;%PL]\n' $input > $output
#HLsamples.tsv文件的形式：CHROM	POS	REF	ALT	VEP_Feature	VEP_SYMBOL	VEP_Gene	GT/sample1 ...

