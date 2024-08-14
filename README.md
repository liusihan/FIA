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

## Citation
If you use BayesQuantify, please cite our paper (thanks!):
> Liu S, Feng X, Bu F. BayesQuantify: an R package utilized to refine the ACMG/AMP criteria according to the Bayesian framework. XXX.


## Getting help
If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussion, please contact Sihan Liu (liusihan@wchscu.cn).
