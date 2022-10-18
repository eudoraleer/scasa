#!/bin/bash
# main parameters
INPUT="/path/to/Test_Dataset"
OUTPUT="/path/to/Test_Dataset/ScasaOut"
ref="/path/to/refMrna.fa.gz"
index="YES" #when index="YES", scasa will index the reference fasta file and write in index_dir. This index_dir cam be reused for other run
index_dir="/path/to/PreBuilt_REF_INDEX" #when index="NO", scasa will use directly the reference indexing in index_dir
nthreads=4
tech="10xv3"
whitelist="/path/to/Test_Dataset/Sample_01_Whitelist.txt"
cellthreshold="none"
project="My_Project"
# other parameters
samplesheet="NULL"
mapper="salmon_alevin"
xmatrix="alevin"
postalign_dir=""
createxmatrix="NO"
