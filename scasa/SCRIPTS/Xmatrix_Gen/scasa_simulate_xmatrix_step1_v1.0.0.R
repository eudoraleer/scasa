#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2023-09-09
##############################################################################
### Nghia/12Sep2023: 
# - speed up using parallel
####

packages <- c("GenomicFeatures","Biostrings")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!all(packages%in% rownames(installed.packages()))){
    needed <- packages[which(!packages %in% rownames(installed.packages()))]
    BiocManager::install(needed)
}

suppressMessages(library("GenomicFeatures"))
suppressMessages(library("Biostrings"))

args = commandArgs(trailingOnly=TRUE)
output_dir <- as.character(args[1])
ref_fasta <- as.character(args[2])
xmatrix_fasta <- as.character(args[3])

library("polyester")
dir.create(output_dir, recursive = T)

tx.all.fasta <- readDNAStringSet(ref_fasta)
tx.export.fasta=tx.all.fasta
writeXStringSet(tx.export.fasta, xmatrix_fasta)

readmat_input <- as.matrix(2*width(tx.export.fasta))
readmat_input[readmat_input<1000] <- 1000

data(cdnaf)
#simulate_experiment_countmat(xmatrix_fasta,
#                             readmat=readmat_input, outdir=output_dir,
#                             readlen = 91, paired = T, seed = 1000,
#                             fraglen = 400,
#                             fragsd = 100,
#                             allow.nonnarrowing = FALSE,
#                             distr = 'normal',
#                             strand_specific = TRUE,
#                             error_rate = 0 , bias="cdnaf")

### simulate by chunks to avoid out of memory

maxNum=length(tx.all.fasta)

chunkSize=10000
chunkNum=trunc((maxNum-1)/chunkSize)+1

r1_fn=paste0(output_dir,"/sample_01_1.fasta")
r2_fn=paste0(output_dir,"/sample_01_2.fasta")

library(foreach)
library(doParallel)
core = detectCores()
if (core > 1) core=round(core/2)
registerDoParallel(cores=core)

res=foreach(chunkInd=1:chunkNum) %dopar% {
#for (chunkInd in 1:chunkNum){
    library("polyester")
    
    chunkStart=(chunkInd-1)*chunkSize+1
    chunkEnd=chunkInd*chunkSize
    if (chunkEnd > maxNum) chunkEnd=maxNum

    xmatrix_fasta_chunk=paste0(xmatrix_fasta,"_chunk_",chunkInd)
    tx.export.fasta_chunk=tx.export.fasta[chunkStart:chunkEnd]
    writeXStringSet(tx.export.fasta_chunk, xmatrix_fasta_chunk)
    readmat_input_chunk <- as.matrix(2*width(tx.export.fasta_chunk))
    readmat_input_chunk[readmat_input_chunk<1000] <- 1000

    output_dir_chunk=paste0(output_dir,"/chunk",chunkInd,"/")

    simulate_experiment_countmat(xmatrix_fasta_chunk,
                             readmat=readmat_input_chunk, outdir=output_dir_chunk,
                             readlen = 91, paired = T, seed = 1000,
                             fraglen = 400,
                             fragsd = 100,
                             allow.nonnarrowing = FALSE,
                             distr = 'normal',
                             strand_specific = TRUE,
                             error_rate = 0 , bias="cdnaf")
    return(NULL)
}

for (chunkInd in 1:chunkNum){
    output_dir_chunk=paste0(output_dir,"/chunk",chunkInd,"/")
    if (chunkInd==1){
        system(paste0("cp ",output_dir_chunk,"/*.fasta ",output_dir))
    }else{
        system(paste0("cat ",output_dir_chunk,"/sample_01_1.fasta >> ",r1_fn))
        system(paste0("cat ",output_dir_chunk,"/sample_01_2.fasta >> ",r2_fn))
    }
    
    system(paste0("rm -r ",output_dir_chunk))
}

system(paste0("rm ",xmatrix_fasta,"_chunk_*"))

xmatrix_dir <- paste(output_dir,"/SIMULATE_XMATRIX/", sep = "")


dir.create(xmatrix_dir, recursive = TRUE)
system(paste("mv ", output_dir,"/sample_01_1.fasta ",xmatrix_dir,"/sample_01_2.fasta", sep = ""))
system(paste("mv ", output_dir,"/sample_01_2.fasta ",xmatrix_dir,"/sample_01_1.fasta", sep = ""))
