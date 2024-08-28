# FIA: An Iterative Analysis Framework for Enhanced Variant Classification in Large Mendelian Disease Cohorts

## Overview
To address these limitations and improve the efficiency and accuracy of genetic diagnosis for Mendelian diseases, we developed the Framework of Iterative Analysis (FIA). This framework integrated the construction of disease mutation databases, automated retrieval of ACMG/AMP evidence, adjustments to specific rules, and the classification of variant pathogenicity. 

Figure 1 illustrates the primary steps involved in FIA.

![Figure 1](https://github.com/liusihan/FIA/blob/main/Figure1.jpeg)
<p align="center"> Figure 1. Schematic overview of FIA. </p>

## Installation

### FIA scripts
Install `FIA` from GitHub:

``` linux
git clone https://github.com/liusihan/FIA.git
```

### Conda
In order to install the dependencies of FIA, we recommend using a dedicated conda environment.
First, you will need the `Conda` Python distribution and package manager. 

```shell
# Download conda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Set permissions to execute
chmod +x Miniconda3-latest-Linux-x86_64.sh 	

# Execute. Make sure to "yes" to add the conda to your PATH
sh ./Miniconda3-latest-Linux-x86_64.sh 		

# Add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

After installing Miniconda, run the following commands to install FIA's dependencies:

```
conda env create -f FIA.yaml
conda activate FIA
pip install -r requirements.txt
```
## Download annotation files
`FIA` use VEP and vcfanno to determine the effect of variants (SNVs, insertions, deletions) on genes, transcripts, and protein sequence. To get annotation for the variant, indexed_vep_cache (homo_sapiens_refseq 105_GRCh37 and 105_GRCh38) and fasta files are required. VEP cache and plugins files can be downloaded as follows:
```shell
# indexed vep cache
mkdir $HOME/.vep
cd $HOME/.vep
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_refseq_vep_105_GRCh38.tar.gz
tar xzf homo_sapiens_refseq_vep_105_GRCh37.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_merged_vep_105_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_merged_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_merged_vep_105_GRCh38.tar.gz
tar xzf homo_sapiens_merged_vep_105_GRCh37.tar.gz
mkdir Plugins
cp -r FIA/Plugins ./

#fasta file
cd FIA/Scripts/AutoPVS1/data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip hg19.fa.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz hg38.fa.gz
gunzip hg38.fa.gz
cat hg19.fa | sed 's/^>chr/>/g' >hg19.nochr.fa
cat hg38.fa | sed 's/^>chr/>/g' >hg38.nochr.fa
samtools faidx hg19.fa
samtools faidx hg38.fa
samtools faidx hg19.nochr.fa
samtools faidx hg38.nochr.fa
```

The database and Scripts directory contain the data and conf for a full example of the GJB2 gene. To run with a complete genome, users needs download the appropriate databases and reverse the file name in the anno.demo.conf file.

## USAGE
```shell
cd FIA
sh FIA_annotation.sh example/1000G_GJB2.vcf.gz GRCh37 Scripts/AutoPVS1/data/hg19.nochr.fa $HOME/.vep test
sh FIA_classification.sh test/annotated_raw.AF_stats.vcf.gz GRCh37 test/FIA
sh test/FIA/bin/run_FIA.sh
```

## Citation
If you use FIA, please cite our paper (thanks!):
> Liu S, et.al,. FIA: An Iterative Analysis Framework for Enhanced Variant Classification in Large Mendelian Disease Cohorts. XXX.


## Getting help
If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussion, please contact Sihan Liu (liusihan@wchscu.cn).
