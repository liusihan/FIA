
samplevcf=$1
genome=$2

scriptdir=$(cd $(dirname ${BASH_SOURCE[0]}); pwd )
#
inputPvariant=$scriptdir/inputfile/PUB_P.txt
inputLinkage=$scriptdir/inputfile/BP2_linkage.txt
genebed=$scriptdir/inputfile/gene_list.txt
pm3_remove=$scriptdir/inputfile/PM3_remove.txt
manual_evidence=$scriptdir/inputfile/Manual_evidence.txt
control_list=$scriptdir/inputfile/control_list.txt
hotregion=$scriptdir/inputfile/PM1_hot_region.bed
hl_samplelist=$scriptdir/inputfile/case.txt
ARgenes=$scriptdir/inputfile/AR_genelist.txt

outdir_name=$3

if [ ! -d "$outdir_name" ]; then
	mkdir $3
fi

outdir=$(realpath $3)

nosampledir_both=$outdir/nonsampledir_both
nosampledir_AR=$outdir/nonsample_evidence_ar
nosampledir_AD=$outdir/nonsample_evidence_ad
nosampledir_other_dir=$outdir/nonsample_evidence_other
sampledir=$outdir/sample_evidence
samplelist=$outdir/samplelist.txt

mkdir -p $outdir/bin $nosampledir_both $nosampledir_AD $nosampledir_AR $sampledir $nosampledir_other_dir 
mkdir -p $outdir/data
mkdir -p $outdir/samplelabling
mkdir -p $outdir/P_satat
mkdir -p $outdir/result
### 脚本路径
combinelabels=$scriptdir/bin/combine_labels_20181207.R
intervar=$scriptdir/bin/intervar_samplelist_20190625_statspecific.pl
statP=$scriptdir/bin/P_allele_stat_20181205.pl

zcat $samplevcf | cut -f1-10 | bgzip >$outdir/data/nonsample.vcf.gz
tabix $outdir/data/nonsample.vcf.gz

nosamplevcf=$outdir/data/nonsample.vcf.gz

echo "steps" > $outdir/bin/record.txt
### 0. data preparing
inputvcf=$(readlink -f $nosamplevcf)
samplevcf=$(readlink -f $samplevcf)
#echo "
#bcftools view -R $genebed $rawvcf -O z > ${samplevcf}.temp.vcf.gz
#bcftools sort ${samplevcf}.temp.vcf.gz -O z > ${samplevcf}
#bcftools index -t ${samplevcf}
#gzip -dc $samplevcf | cut -f1-9 > $outdir/data/ACMG_selectedgene_nosample.vcf
#bcftools sort $outdir/data/ACMG_selectedgene_nosample.vcf -T /public/home/swgenetics_3/zengyuanyuan/temp/|bgzip -c > $inputvcf
#tabix -f -p vcf $inputvcf
#" > $outdir/bin/data_preparing.00.sh

###  1. evidence labling

## 
echo "
python $scriptdir/Scripts/FIA_PVS1.py -i $inputvcf -g $genome -o ${nosampledir_both}/
perl $scriptdir/bin/addPS1_PM5_by_transcript_20181011.pl -anno $inputvcf -t $inputPvariant -o ${nosampledir_both}/PS1_PM5.tsv
perl $scriptdir/bin/addBP2_by_transcript_20181011.pl -anno $inputvcf -t $inputPvariant -l $inputLinkage -o ${nosampledir_both}/BP2.tsv
##PP3,BP4
sh $scriptdir/bin/pp3_bp4.sh $inputvcf $outdir

##
sh $scriptdir/bin/execute_v4_change_pp3_bp4.sh $inputvcf $inputPvariant ${nosampledir_AD} ${nosampledir_AR}
##PM1
sh $scriptdir/bin/pm1.sh $inputvcf ${nosampledir_other_dir} ${hotregion}

echo \"step1.1 finished\" >> $outdir/bin/record.txt
" > $outdir/bin/generate_evidence.0101.sh
##
echo "
sh $scriptdir/bin/bs2.sh $samplevcf ${nosampledir_other_dir} $inputPvariant ${control_list}
echo \"step1.2 finished\" >>$outdir/bin/record.txt
" > $outdir/bin/generate_evidence.0102.sh
## PM3
echo "sh $scriptdir/bin/pm3.sh $samplevcf ${nosampledir_other_dir} $inputPvariant $pm3_remove $hl_samplelist $ARgenes
echo \"step1.3 finished\" >> $outdir/bin/record.txt
" > $outdir/bin/generate_evidence.0103.sh
## PP4 
echo "
mkdir -p $sampledir/rawData/
mkdir -p $sampledir/pp4Result/
cp $scriptdir/inputfile/geneset-*.txt $sampledir/rawData/
cp $scriptdir/inputfile/phenotype.txt $sampledir/rawData/

sh $scriptdir/bin/pp4/getHL_122genes.sh $sampledir/rawData/ $samplevcf
perl $scriptdir/bin/pp4/reformat.pl $sampledir/rawData/ HLsamples_122genes.tsv HLsamples_122genes.reformat.tsv
python3.6 $scriptdir/bin/pp4/pp4_main.py $sampledir/rawData/phenotype.txt $sampledir/rawData/geneset-*.txt \"$sampledir/pp4Result/\" \"$sampledir/rawData/HLsamples_122genes.reformat.tsv\" \"$samplevcf\"
echo \"step1.4 finished\" >> $outdir/bin/record.txt
" > $outdir/bin/generate_evidence.0104.sh




### 2. combine evidence
echo "
paste $nosampledir_both/*.tsv | cut -f1-6,12-13,19 > $outdir/none_sample_evidence_both.tsv
paste $nosampledir_AD/*.tsv > $outdir/none_sample_evidence_AD.tsv
paste $nosampledir_AR/*.tsv > $outdir/none_sample_evidence_AR.tsv
echo \"step2 finished\" >> $outdir/bin/record.txt
" >$outdir/bin/combine_nonsample_evidence.02.sh

### 3. intervar

echo "cp $scriptdir/inputfile/all_prediction_20181105.txt $outdir/
R --vanilla --no-echo -f $combinelabels --args $outdir/none_sample_evidence_both.tsv $outdir/none_sample_evidence_AR.tsv $outdir/none_sample_evidence_AD.tsv ${nosampledir_other_dir}
ls $sampledir/pp4Result/*.tsv > $samplelist
#mv $outdir/ /NAS/shliu/
#ln -s /NAS/shliu/$outdir_name $outdir 
perl $intervar $outdir/data/test_dat_ar.txt $outdir/data/test_dat_ad.txt $samplelist ${manual_evidence} $outdir/samplelabling/
## variants filtering 
sh $scriptdir/bin/ACMG_filter_sites.sh $inputvcf $manual_evidence $outdir
### stat
#Rscript $scriptdir/bin/evidence_split.R $outdir/result/ $outdir/result/ACMG_classfy_results_filtered.forstat.tsv $outdir/result/ACMG_split_feature.tsv
#mkdir -p $outdir/result/stat/
#Rscript $scriptdir/bin/stat.R $outdir/result/stat/ $outdir/result/ /public/home/swgenetics_4/database/InputFiles/P_LP_SNP.txt
" >$outdir/bin/intervar_samplelist.sh


###### 4. submit to SC

echo "
nohup sh $outdir/bin/generate_evidence.0101.sh > $outdir/bin/generate_evidence.0101.log 2>&1 &
nohup sh $outdir/bin/generate_evidence.0102.sh > $outdir/bin/generate_evidence.0102.log 2>&1 &
nohup sh $outdir/bin/generate_evidence.0103.sh > $outdir/bin/generate_evidence.0103.log 2>&1 &
nohup sh $outdir/bin/generate_evidence.0104.sh > $outdir/bin/generate_evidence.0104.log 2>&1 &

for i in {1..900}
do
	record=\`awk 'END{print NR}' $outdir/bin/record.txt\`
	if [ \"\$record\" -eq 5 ]
	then
		break
	else
		sleep 2m
	fi
done

nohup sh $outdir/bin/combine_nonsample_evidence.02.sh > $outdir/bin/combine_nonsample_evidence.02.log 2>&1 &

for i in {1..1000}
do
    record=\`awk 'END{print NR}' $outdir/bin/record.txt\`
    if [ \"\$record\" -eq 6 ]
    then
        break
    else
        sleep 30s
    fi
done

nohup sh $outdir/bin/intervar_samplelist.sh > $outdir/bin/intervar_samplelist.log 2>&1 &

" > $outdir/bin/run_FIA.sh
