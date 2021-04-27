<img alt="scasa logo" src="https://github.com/eudoraleer/scasa/blob/main/doc/SCASA_LOGO.png">

# scasa
Single cell transcript quantification tool

__scasa__ is a single cell transcript quantification software designed for single cell RNA-Sequencing data. The software comprises pseudo-alignment to quantification steps. See the [__scasa &#124; wiki__](https://github.com/eudoraleer/scasa/wiki) for more details on __scasa__.

## Installation

Scasa only has Linux version at the moment. The software is already compiled for Linux and installation time needed is less than two minutes. Please follow the instructions below to install scasa:

Dependency packages below are needed to run the quick tutorial below for scasa:

##### (1) [__Alevin__](https://salmon.readthedocs.io/en/latest/alevin.html)

##### (2) [__R__](https://www.r-project.org)

After dependency packages are installed:

1. Download scasa from our Github scasa release  [__scasa__](https://github.com/eudoraleer/scasa/releases/tag/scasa.v1.0.0) and untar the downloaded file:
```sh
wget https://github.com/eudoraleer/scasa/releases/download/scasa.v1.0.0/scasa_v1.0.0.tar.gz
tar -xzvf scasa_v1.0.0.tar.gz
```

2. Add __scasa__ folder to environment variables PATH:
```sh
export PATH=$PWD/scasa:$PATH
```

3. Now you are ready to use scasa!

## Quick tutorial on scasa

After installation, test out scasa by typing  `scasa --help`  in the terminal to see a list of available commands. To see a list of detailed options on scasa, visit our [__wiki__](https://github.com/eudoraleer/scasa/wiki) page.

1. Download our [__Test_Dataset__](https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz) and unzip it:
```sh
wget https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz
tar xvzf Test_Dataset.tar.gz              
```

2. Download the cDNA fasta of [__hg38: refMrna__](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/)
```sh
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz
```

3. Enter the following command to kick start the analysis (set a higher number of threads to enable faster processing):
```sh
cd Test_Dataset
scasa --fastq Sample_01_S1_L001_R1_001.fastq,Sample_01_S1_L001_R2_001.fastq \
      --ref <hg38_ref_file_path>  \
      --whitelist <test_dataset_whitelist_path> \
      --nthreads 4
```
4. After completed your analysis with scasa, you will see that scasa has generated a project output directory with name `<SCASA_project_name_timestamp>` with the following sub-directories:

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

### A Quick Example:
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
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz
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
### Example of isoform quantification with Scasa for the CITE-seq RNA samples of [__Bone Marrow Mononuclear Cells__](https://www.ebi.ac.uk/ena/browser/view/PRJNA528319):
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
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz
refPath=$PWD/refMrna.fa.gz
cd ..

##################################################################
# 4. Download test dataset:
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

