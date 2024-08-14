#cd /public/home/swgenetics_3/zengyuanyuan/backup/ACMG_zmj
input=$1
out_dir=$2
hotregion=$3
#cat header.tsv >${out_dir}/PM1.tsv
echo -e "chrom\tpos\tref\talt\ttranscript\tPM1" >${out_dir}/PM1.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VEP_Feature\t2\n' -R $hotregion $input | perl -F'\t' -alne '@fields=split(/\|/,$F[4]);for $each (@fields) { print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$each\t$F[5]";}' >> ${out_dir}/PM1.tsv
