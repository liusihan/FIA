# FIA: An iterative analysis framework for semi-automating variant classification in genetic hearing loss

## Overview
Here, we propose a 'framework of iterative analysis (FIA)' designed to semi-automate variant classification in genetic hearing loss, particularly suited for cohort data. The development of FIA aims to improve the efficiency and accuracy of genetic diagnosis across various deafness genes, with potential applicability to other Mendelian disorders. FIA crucially leverages genetic and phenotypic connections within the patient cohort, utilizing additional evidence from diagnosed cases to iteratively identify novel P/LP variants in previously undiagnosed cases. 

The main steps of FIA are shown in Figure 1:

## Installation

### FIA scripts
Install `FIA` from GitHub:

``` linux
git clone https://github.com/liusihan/FIA
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
## Annotation
`FIA` use VEP and vcfanno to determine the effect of variants (SNVs, insertions, deletions) on genes, transcripts, and protein sequence. To get annotation for the variant, indexed_vep_cache (homo_sapiens_refseq 105_GRCh37 and 105_GRCh38) and fasta files are required. VEP cache and faste files can be downloaded as follows:
```shell
# indexed vep cache
cd $HOME/.vep
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_refseq_vep_105_GRCh38.tar.gz
tar xzf homo_sapiens_refseq_vep_105_GRCh37.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_merged_vep_105_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-105/variation/vep/homo_sapiens_merged_vep_105_GRCh37.tar.gz
tar xzf homo_sapiens_merged_vep_105_GRCh38.tar.gz
tar xzf homo_sapiens_merged_vep_105_GRCh37.tar.gz

#fasta file
cd FIA/Scripts/AutoPVS1/data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip hg19.fa.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg19.fa
samtools faidx hg38.fa
```

## Citation
If you use BayesQuantify, please cite our paper (thanks!):
> Liu S, Feng X, Bu F. BayesQuantify: an R package utilized to refine the ACMG/AMP criteria according to the Bayesian framework. XXX.


## Getting help
If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussion, please contact Sihan Liu (liusihan@wchscu.cn).
