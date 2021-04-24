<img alt="scasa logo" src="https://github.com/eudoraleer/scasa/blob/main/doc/SCASA_LOGO.png">

# scasa
Single cell transcript quantification tool

__scasa__ is a single cell transcript quantification software designed for single cell RNA-Sequencing data. The software comprises pseudo-alignment to quantification steps. See the [__scasa &#124; wiki__](https://github.com/eudoraleer/scasa/wiki) for detailed examples and instructions on how to use __scasa__ as part of a single-cell RNA-seq workflow.

## Installation

Scasa only has Linux version at the moment. To install, download the scasa folder from the. The software is already compiled for Linux and installation time needed is less than two minutes.

Dependency packages below are needed to run the quick tutorial below for scasa:

##### (1) [__Alevin__](https://salmon.readthedocs.io/en/latest/alevin.html)

##### (2) [__R__](https://www.r-project.org)

##### (3) [__BBMap__](https://github.com/BioInfoTools/BBMap)

After depency packages are installed:

1. Download scasa folder from our Github repository  [__scasa__](https://github.com/eudoraleer/scasa) and do the following:

        git clone https://github.com/scasa/

2. Add __scasa__ folder to environment variables PATH:

        PATH=<scasa directory>:$PATH
        
3. Now you are ready to use scasa!

## Quick tutorial on scasa

After installation, test out scasa by typing  `scasa --help`  in the terminal to see a list of available commands. To see a list of detailed options on scasa, visit our [__wiki__](https://github.com/eudoraleer/scasa/wiki) page.

1. Download our provided [__Test_Dataset__](https://github.com/eudoraleer/scasa/tree/main/Test_Dataset) from our Github repository (if you have downloaded you can skip this step).

2. Download the following files:

    [__Hg38 reference: hg38 refMrna__](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/refMrna.fa.gz)

2. Enter the following command to kick start the analysis (set a higher number of threads to enable faster processing):

        cd <Test_Dataset_directory>
        
        scasa --fastq Sample_01_S1_L001_R1_001.fastq,Sample_01_S1_L001_R2_001.fastq \
              --ref <hg38_ref_file_path>  \
              --whitelist <test_dataset_whitelist_path> \
              --nthreads 4


