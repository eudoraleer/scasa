<img alt="scasa logo" src="https://github.com/eudoraleer/scasa/blob/main/doc/SCASA_LOGO.png">

# Scasa

Single cell transcript quantification tool

__Scasa__ is a single cell transcript quantification tool tailored for single cell RNA-Sequencing data. The software comprises pseudo-alignment to quantification steps. See the [__scasa &#124; wiki__](https://github.com/eudoraleer/scasa/wiki) for more details on __scasa__.

##### (Datasets and codes to reproduce the results in the __scasa__ study could be view from Scasa Paper GitHub Page: [__Scasa_Paper__](https://github.com/eudoraleer/Scasa_Paper))

If you are using Scasa in your research, please cite:

Lu Pan, Huy Q Dinh, Yudi Pawitan, Trung Nghia Vu, Isoform-level quantification for single-cell RNA sequencing, Bioinformatics, 2021;, btab807, https://doi.org/10.1093/bioinformatics/btab807

## Installation

Scasa only has Linux version at the moment. The software is already compiled for Linux and installation time needed is less than two minutes. Please follow the instructions below to install scasa:

Dependency packages below are needed to run the quick tutorial below for scasa:

##### (1) [__Alevin__](https://salmon.readthedocs.io/en/latest/alevin.html)

##### (2) [__R__](https://www.r-project.org)
```sh
# R packages needed (they will be automatically installed by scasa provided that enough permission is given to install R packages under default R library directory):
library(GenomicFeatures)
library(Biostrings)
library(polyester)
library(foreach)
library(doParallel)
library(data.table)
library(plyr)
```

After dependency packages are installed:

1. Download scasa from our Github scasa release  [__scasa__](https://github.com/eudoraleer/scasa/releases/tag/scasa.v1.0.0) and untar the downloaded file:
```sh
wget https://github.com/eudoraleer/scasa/releases/download/scasa.v1.0.1/scasa_v1.0.1.tar.gz
tar -xzvf scasa_v1.0.1.tar.gz
```

2. Add __scasa__ folder to environment variables PATH:
```sh
export PATH=$PWD/scasa:$PATH
```

3. Now you are ready to use scasa!

## Quick Tutorial on Scasa

After installation, test out scasa by typing  `scasa --help`  in the terminal to see a list of available commands. To see a list of detailed options on scasa, visit our [__wiki__](https://github.com/eudoraleer/scasa/wiki) page.

1. Download our [__Test_Dataset (200 cells)__](https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz) and unzip it:
```sh
wget https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz
tar xvzf Test_Dataset.tar.gz              
```

2. Download the cDNA fasta of [__hg38: refMrna__](https://www.dropbox.com/s/xoa6yl562a5lv35/refMrna.fa.gz)
#### NOTE: The current pre-built X-maxtrix was built based on this refMrna.fa.gz file (downloaded from ucsc goldenPath - 16Feb2021). There might have slight differences between this pre-built version and other version of hg38-refMrna which can break the running. 
#### If you want to other reference annotations, please follow the instruction in the [__wiki__](https://github.com/eudoraleer/scasa/wiki) page and [How-to-run-Scasa-for-a-new-annotation](https://github.com/eudoraleer/scasa/wiki/How-to-run-Scasa-for-a-new-annotation) to rebuid then replace the X-matrix and other files in scasa/REFERENCE folder.

3. Enter the following command to kick start the analysis (set a higher number of threads to enable faster processing):
```sh
cd Test_Dataset
scasa --fastq Sample_01_S1_L001_R1_001.fastq,Sample_01_S1_L001_R2_001.fastq \
      --ref <hg38_ref_file_path>  \
      --whitelist <test_dataset_whitelist_path> \
      --nthreads 4
```
4. After you have completed your analysis with scasa, you will see that scasa has generated a project output directory with name `<SCASA_project_name_timestamp>` with the following sub-directories:

        <SCASA_project_name_timestamp>/
        ├── LOG/
        ├── 0PRESETS/
        ├── 1ALIGN/
        └── 2QUANT/
            ├──<sample_1_quantification_output>
            │   └──scasa_isoform_expression.txt
            │   └──scasa_gene_expression.txt
            └──..

5. Isoform and gene expression output can be found under the `2QUANT/` directory in the output folder:
```sh
cd <SCASA_project_name_timestamp>/2QUANT/<sample_1_quantification_output>/
```
Now that you have learnt how to run scasa!

## Annotation
In the scasa study, we used the Refseq hg38 for all analysis. This annotation does not contain the transcripts of mitochondria chromosome which are sometimes required for some analyses. Therefore we added the transcripts of the chM from GRCh38 annotation to build a new reference data for scasa.

The annotation can be downloaded here: [__Annotation data for running scasa (alevin mapper) using refseq hg38 with chM__](https://kise-my.sharepoint.com/:u:/r/personal/trungnghia_vu_ki_se/Documents/Shared_folder/Public/scasa_refseq_hg38_chrM.zip?csf=1&web=1&e=b2Iebj). The version for Homo_sapiens.GRCh38.106 of ENSEMBL can be found in [folder Anno](https://github.com/eudoraleer/scasa/tree/main/Anno).
Users can use the annotation files and replace the default reference data in scasa to get results with chM.

It is noted that we have not tested the performances of scasa for other annotation systems such as ENCODE/ENSEMBL.

For those who want to run scasa for new annotations, some scripts in [__aux folder__](https://github.com/eudoraleer/scasa/tree/main/aux) are helpful to generate Xmatrix and necessary reference data. See more instructions in [How-to-run-Scasa-for-a-new-annotation](https://github.com/eudoraleer/scasa/wiki/How-to-run-Scasa-for-a-new-annotation).

## Scasa Quantification Run Time (Hours)

On Linux CentOs 7, we tested from thread number 1 to 64 for both [__small simulated dataset__](https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz) (__200 cells, dataset from Step 1, Quick tutorial on scasa__) and for a __larger simulated dataset (3955 cells)__ and below are the runtime information for both simulated datasets in terms of hours:

<img alt="scasa runtime" src="https://github.com/eudoraleer/scasa/blob/main/doc/SCASA_QUANTIFICATION_RUN_TIME_HOURS.png">


## A Quick Example on The Small Simulated Dataset (200 Cells)
### (Dataset from __Step 1, Quick Tutorial on Scasa__)
```sh
##################################################################
# 1. Download scasa:
##################################################################
wget https://github.com/eudoraleer/scasa/releases/download/scasa.v1.0.0/scasa_v1.0.0.tar.gz
tar -xzvf scasa_v1.0.0.tar.gz
export PATH=$PWD/scasa:$PATH

##################################################################
# 2. Download salmon alevin:
##################################################################
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz
tar -xzvf salmon-1.4.0_linux_x86_64.tar.gz
export PATH=$PWD/salmon-latest_linux_x86_64/bin:$PATH
export LD_LIBRARY_PATH=$PWD/salmon-latest_linux_x86_64/lib:$LD_LIBRARY_PATH

##################################################################
# 3. Download UCSC hg38 cDNA fasta reference:
##################################################################
mkdir Annotation
cd Annotation
wget https://www.dropbox.com/s/xoa6yl562a5lv35/refMrna.fa.gz
refPath=$PWD/refMrna.fa.gz
cd ..

##################################################################
# 4. Download test dataset:
##################################################################
wget https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz
tar xvzf Test_Dataset.tar.gz
cd Test_Dataset

##################################################################
# 5. Run scasa:
##################################################################
scasa --fastq Sample_01_S1_L001_R1_001.fastq,Sample_01_S1_L001_R2_001.fastq \
      --ref $refPath \
      --whitelist Sample_01_Whitelist.txt \
      --nthreads 2 \
      --out Scasa_out
      
##################################################################
# DONE!
##################################################################

```

## Using docker to run scasa
Using docker can avoid the issues of the installation of the dependent tools and the enviroment. Users need to revise the docker_params.sh to the paths to input and ouput folders and other scasa parameter setting. The scripts below show an example of running docker for the test data above.
```sh
#1) Pull the docker of scasa to use:
sudo docker pull nghiavtr/scasa:v1.0.1

#2) Download test dataset:
wget https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz
tar xvzf Test_Dataset.tar.gz

#3) Download UCSC hg38 cDNA fasta reference:
wget https://www.dropbox.com/s/xoa6yl562a5lv35/refMrna.fa.gz

#4) Download runScasaDocker.sh:
wget https://raw.githubusercontent.com/eudoraleer/scasa/main/docker/runScasaDocker.sh

#5) Download docker_params.sh:
wget https://raw.githubusercontent.com/eudoraleer/scasa/main/docker/docker_params.sh

#6) replace "/path/to/" in the docker_params.sh by the current path
# Please revise the paths in the docker_params.sh following your project
istr="/path/to"
ostr="$PWD"
eval "sed -i -e 's#"$istr"#"$ostr"#g' docker_params.sh"

#7) Run scasa using docker with the parameter settings in docker_params.sh
sudo bash runScasaDocker.sh -param docker_params.sh
```
## Real Sample Example on [__Bone Marrow Mononuclear Cells__](https://www.ebi.ac.uk/ena/browser/view/PRJNA528319) CITE-seq Dataset
```sh
##################################################################
# 1. Download scasa:
##################################################################
wget https://github.com/eudoraleer/scasa/releases/download/scasa.v1.0.0/scasa_v1.0.0.tar.gz
tar -xzvf scasa_v1.0.0.tar.gz
export PATH=$PWD/scasa:$PATH

##################################################################
# 2. Download salmon alevin:
##################################################################
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz
tar -xzvf salmon-1.4.0_linux_x86_64.tar.gz
export PATH=$PWD/salmon-latest_linux_x86_64/bin:$PATH
export LD_LIBRARY_PATH=$PWD/salmon-latest_linux_x86_64/lib:$LD_LIBRARY_PATH

##################################################################
# 3. Download UCSC hg38 cDNA fasta reference:
##################################################################
mkdir Annotation
cd Annotation
wget https://www.dropbox.com/s/xoa6yl562a5lv35/refMrna.fa.gz
refPath=$PWD/refMrna.fa.gz
cd ..

##################################################################
# 4. Download the CITE-seq RNA samples:
##################################################################

mkdir CiteSeqData
InputDir=$PWD/CiteSeqData
cd CiteSeqData
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/003/SRR8758323/SRR8758323_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/003/SRR8758323/SRR8758323_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/004/SRR8758324/SRR8758324_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/004/SRR8758324/SRR8758324_2.fastq.gz
cat SRR8758323_1.fastq.gz SRR8758324_1.fastq.gz > HBMC_Stuart2019_RNA_L001_R1_001.fastq.gz
cat SRR8758323_2.fastq.gz SRR8758324_2.fastq.gz > HBMC_Stuart2019_RNA_L001_R2_001.fastq.gz
rm SRR8758323_1.fastq.gz SRR8758323_2.fastq.gz SRR8758324_1.fastq.gz SRR8758324_2.fastq.gz

##################################################################
# 5. Run scasa:
##################################################################

threadNum=16
scasa --in $InputDir --fastq HBMC_Stuart2019_RNA_L001_R1_001.fastq.gz,HBMC_Stuart2019_RNA_L001_R2_001.fastq.gz --ref $refPath --cellthreshold 35000 --tech 10xv2 --nthreads $threadNum --out Scasa_out

#################################################################
# DONE!
##################################################################

```

