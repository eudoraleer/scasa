#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 3_1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
packages <- c("data.table")

if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages, repos = "https://cran.r-project.org")
}

library("data.table")

args = commandArgs(trailingOnly=TRUE)
input_dir <- as.character(args[1])
TxCellID_mapping_fn <- as.character(args[2])
output_dir <- as.character(args[3])

bus_output <- fread(paste(input_dir,"/output.bus.txt", sep = ""))
colnames(bus_output) <- c("Barcode","UMI_Code","eqClass","Count")
load(TxCellID_mapping_fn)

eq <- read.table(paste(input_dir,"/matrix.ec", sep = ""),stringsAsFactors=FALSE)
colnames(eq) <- c("eqClass","Transcript_Number")
eq$eqClass <- as.character(eq$eqClass)
eq$Transcript_Number <- as.character(eq$Transcript_Number)

ref <- read.table(paste(input_dir,"/transcripts.txt", sep = ""))
colnames(ref) <-c("Transcript_ID")
ref$Transcript_Number <- seq(0,(nrow(ref)-1),1)
#25Jan2021/Nghia: fix the bug when a total number of transcripts is greater than 100000(1e+05)
ref$Transcript_Number = as.character(as.integer(ref$Transcript_Number))

s <- strsplit(eq$Transcript_Number, split = ",")
eq <- data.frame(eqClass = rep(eq$eqClass, sapply(s, length)), Transcript_Number = unlist(s))
eq$Transcript_ID <- ref[match(eq$Transcript_Number,ref$Transcript_Number),"Transcript_ID"]

#### get count
mappedReads=bus_output[, .N, by = eqClass]

eq$eqClass=as.integer(as.character(eq$eqClass))
allEq=unique(eq$eqClass)
pick=mappedReads$eqClass %in% allEq
sum(pick)
dim(mappedReads)
  
pick=which(eq$eqClass %in% mappedReads$eqClass)
current_final=eq[pick,c("Transcript_ID","eqClass")]
length(unique(current_final$eqClass))
current_final=as.data.frame(current_final,stringsAsFactors=FALSE)
current_final$Count=mappedReads$N[match(current_final$eqClass,mappedReads$eqClass)]

#### now we need to get the weight
cell_eq <- fread(paste(input_dir,"/output.barcode_eq.bus.txt", sep = ""))
colnames(cell_eq)="barcode_eq"
txweight=cell_eq[, .N, by = barcode_eq]
#cell_eq=paste(bus_output$Barcode, bus_output$eqClass, sep = "")
#txweight=table(cell_eq)
cellID.txweight=sapply(txweight$barcode_eq,function(x) substr(x,1,16))
names(cellID.txweight)=NULL
matchID=match(cellID.txweight,white_list.txID.all)
TxID.txweight=txID.all[matchID]

eq.txweight=sapply(txweight$barcode_eq,function(x) substring(x,17))
names(eq.txweight)=NULL
eq.txweight=as.integer(eq.txweight)

tx_eq.txweight=paste0(TxID.txweight,eq.txweight)

###update weight to the eqclass table
tx_eq.current_final=paste0(current_final$Transcript_ID,current_final$eqClass)
current_final$Weight <- 0
matchID=match(tx_eq.txweight,tx_eq.current_final)
#25Jan2021/Nghia: fix the bug when there are new combinations of (tx,eqc)
pick=is.na(matchID)
sum(txweight$N[pick])
current_final$Weight[matchID[!pick]]=txweight$N[!pick]
write.table(current_final,paste(output_dir,"/Xmatrix_eqClass.txt", sep = ""), quote = F, row.names = F,sep="\t")
