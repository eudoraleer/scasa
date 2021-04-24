#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
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
simulate_experiment_countmat(xmatrix_fasta,
                             readmat=readmat_input, outdir=output_dir,
                             readlen = 91, paired = T, seed = 1000,
                             fraglen = 400,
                             fragsd = 100,
                             allow.nonnarrowing = FALSE,
                             distr = 'normal',
                             strand_specific = TRUE,
                             error_rate = 0 , bias="cdnaf")

xmatrix_dir <- paste(output_dir,"SIMULATE_XMATRIX/", sep = "")
dir.create(xmatrix_dir, recursive = TRUE)
system(paste("mv ", output_dir,"/sample_01_1.fasta ",xmatrix_dir,"/sample_01_2.fasta", sep = ""))
system(paste("mv ", output_dir,"/sample_01_2.fasta ",xmatrix_dir,"/sample_01_1.fasta", sep = ""))
