<img alt="scasa logo" src="https://github.com/eudoraleer/scasa/blob/main/doc/SCASA_LOGO.png">

# scasa
Single cell transcript quantification tool

__scasa__ is a single cell transcript quantification software designed for single cell RNA-Sequencing data. The software comprises pseudo-alignment to quantification steps. See the [__scasa &#124; wiki__](https://github.com/eudoraleer/scasa/wiki) for detailed examples and instructions on how to use __scasa__ as part of a single-cell RNA-seq workflow.

## Installation

Scasa only has Linux version at the moment. To install, download the scasa folder from the. The software is already compiled for Linux and installation time needed is less than two minutes.

Dependency packages below are needed to run the quick tutorial below for scasa:

##### (1) [__Alevin__](https://salmon.readthedocs.io/en/latest/alevin.html)

##### (2) [__R__](https://www.r-project.org) with following installed packages: data.table, foreach, doParallel, GenomicFeatures, Biostrings (scasa will automatically install these R packages while if users have not installed them, but users could also install in prior:)


After dependency packages are installed:

1. Download scasa folder from our Github repository  [__scasa__](https://github.com/eudoraleer/scasa) and do the following:
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

2. Download the cDNA fasta of [__hg38: refMrna__](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz)
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


### Example of scripts to install and run scasa for the test dataset

```sh
#download scasa
wget https://github.com/eudoraleer/scasa/releases/download/scasa.v1.0.0/scasa_v1.0.0.tar.gz
tar -xzvf scasa_v1.0.0.tar.gz
export PATH=$PWD/scasa:$PATH

#download salmon/alevin
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz
tar -xzvf salmon-1.4.0_linux_x86_64.tar.gz
export PATH=$PWD/salmon-latest_linux_x86_64/bin:$PATH
export LD_LIBRARY_PATH=$PWD/salmon-latest_linux_x86_64/lib:$LD_LIBRARY_PATH

#download test samples
wget https://www.dropbox.com/s/gsi8x4fshbn0p11/Test_Dataset.tar.gz
tar xvzf Test_Dataset.tar.gz
cd Test_Dataset

#download cDNA fasta of hg38
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz

#run scasa
scasa --fastq Sample_01_S1_L001_R1_001.fastq,Sample_01_S1_L001_R2_001.fastq --ref refMrna.fa.gz --whitelist Sample_01_Whitelist.txt --nthreads 4 --out Scasa_out

#done

```
