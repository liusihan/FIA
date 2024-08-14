#!/usr/bin/bash

##test for hg19+GRCh37 version（without 'chr'，线粒体写作M）!!!

#####tsv to cf
# SOURCE="input"
# ./tab2vcf --chr 1 \
#          --pos 2 \
#          --ref 3 \
#          --alt 4 \
#          --source $SOURCE \
#          --reference hg38 \
#          --prefix SWHG \
#          validate_sample.tsv \
#          > input.vcf

RawVCF=$(readlink -f $1)
genome=$2
Ref=$(readlink -f $3)
Cache=$(readlink -f $4)
OutputDIR=$(readlink -f $5)

mkdir -p $OutputDIR

chmod -R 777 $OutputDIR

paral_num=10

#####参考基因组版本选择
reference=$Ref

##Cache
CacheDIR=$Cache

##ScriptDIR
ScriptDIR=$(dirname $(readlink -f "$0"))/Scripts

##Inputfile
InputVCF=$RawVCF

FIAdir=$(dirname $(readlink -f "$0"))

cd $OutputDIR

####################Split and Normalize####################

chmod 777 $InputVCF

gunzip -f $InputVCF

InputVCF=${InputVCF%.*}

sed -i 's/AD,Number=./AD,Number=R/' $InputVCF

bgzip -f -@ $paral_num $InputVCF

InputVCF=$InputVCF.gz

bcftools index --threads $paral_num -f -t $InputVCF

bcftools norm --threads $paral_num -f $reference -c w -m -both $InputVCF -Ov > $InputVCF.norm

bcftools sort $InputVCF.norm |sed -e 's/^M/MT/' -e 's/ID=M,/ID=MT,/' > $InputVCF.norm.sort

bgzip -f -@ $paral_num $InputVCF.norm.sort

bcftools index --threads $paral_num -f -t $InputVCF.norm.sort.gz

rm -rf $InputVCF.norm


####################VEP annotation####################
InputVCF=$InputVCF.norm.sort.gz

chmod -R 777 $OutputDIR

chmod 777 $InputVCF $InputVCF.*

time vep \
-cache \
--dir $CacheDIR \
--cache_version 105 \
--offline \
--assembly $genome \
--format vcf \
--fork $paral_num \
--fasta $reference \
-i $InputVCF \
-o $OutputDIR/vep_output.vcf \
--merged \
--force_overwrite \
--everything \
--flag_pick \
--vcf \
--plugin Downstream \
--plugin GeneSplicer,$CacheDIR/Plugins/GeneSplicer/sources/genesplicer,$CacheDIR/Plugins/GeneSplicer/human,tmpdir=$CacheDIR/Plugins/GeneSplicer/mytmp \
--plugin MaxEntScan,$CacheDIR/Plugins/MaxEntScan/fordownload \
--plugin SpliceRegion,Extended



####################Split csq####################
time python $ScriptDIR/split_csq.py $OutputDIR/vep_output.vcf > $OutputDIR/vep_split_csq.MT.vcf
##rm vep_output.vcf

sed -e 's/^MT/M/' -e 's/ID=MT,/ID=M,/' $OutputDIR/vep_split_csq.MT.vcf > $OutputDIR/vep_split_csq.vcf
rm -rf $OutputDIR/vep_split_csq.MT.vcf

####################vcfanno####################
sed -e 's|\~/FIA_github|'$FIAdir'|g' $ScriptDIR/anno.demo.conf > $ScriptDIR/anno.conf

time vcfanno -p $paral_num \
 -lua $ScriptDIR/custom.lua \
 $ScriptDIR/anno.conf \
 $OutputDIR/vep_split_csq.vcf \
> $OutputDIR/annotated_raw.AF_stats.vcf

time bgzip -f -@ $paral_num $OutputDIR/annotated_raw.AF_stats.vcf

time bcftools index --threads $paral_num -f -t $OutputDIR/annotated_raw.AF_stats.vcf.gz

##rm vep_split_csq.vcf
