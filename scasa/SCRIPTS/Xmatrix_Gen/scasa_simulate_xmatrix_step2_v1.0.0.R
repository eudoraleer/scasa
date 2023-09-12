#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 2
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2023-09-09
##############################################################################
### Nghia/12Sep2023: 
# - speed up using parallel
# - optimize memory usage
####


library(data.table)
args = commandArgs(trailingOnly=TRUE)
whitelist <- as.character(args[1])
output_dir <- as.character(args[2])
odd_sim_fasta <- as.character(args[3])

#white_list <- read.table(whitelist, header = F,stringsAsFactors=FALSE)
#white_list <- white_list[,1]
white_list <- fread(whitelist, header = F)
white_list <- white_list$V1


#txread_map=read.table(odd_sim_fasta, header = F,stringsAsFactors=FALSE)
#txread_map=txread_map[,1]
txread_map=fread(odd_sim_fasta, header = F)
txread_map=txread_map$V1

txID.all=unique(txread_map)
num_reads=length(txread_map)

set.seed(2021)
white_list.txID.all=sample(white_list,length(txID.all),replace = FALSE)
white_list=NULL

save(txID.all,white_list.txID.all, file=paste(output_dir,"/TxCellID_mapping.RData",sep = ""))



matchID=match(txread_map,txID.all)
txread_map=NULL
txID.all=NULL

current_white_list=white_list.txID.all[matchID]
write.table(current_white_list, paste(output_dir,"/sample_01.white.list.txt",sep = ""), quote = F, row.names = F, col.names = F)
white_list.txID.all=NULL
current_white_list=NULL
matchID=NULL


print("Finished generating sample white list..")


# UMI,12bp
umi_size <- 12
set.seed(2021)
#current_umi <- do.call(paste0, replicate(umi_size, sample(c("A","T","G","C"), num_reads, TRUE), FALSE))
#write.table(current_umi, paste(output_dir,"/sample_01.umi.txt",sep = ""), quote = F, row.names = F, col.names = F)

#do in chunks to reduce the memory
maxNum=num_reads
chunkSize=1000000
chunkNum=trunc((maxNum-1)/chunkSize)+1

library(foreach)
library(doParallel)
core = detectCores()
if (core > 1) core=round(core/2)
registerDoParallel(cores=core)

res=foreach(chunkInd=1:chunkNum) %dopar% {
#for (chunkInd in 1:chunkNum){
	#cat(" ",chunkInd)
    chunkStart=(chunkInd-1)*chunkSize+1
    chunkEnd=chunkInd*chunkSize
    if (chunkEnd > maxNum) chunkEnd=maxNum

    num_reads_chunk=chunkEnd-chunkStart+1

	current_umi <- do.call(paste0, replicate(umi_size, sample(c("A","T","G","C"), num_reads_chunk, TRUE), FALSE))
	write.table(current_umi, paste0(output_dir,"/sample_01.umi_chunk_",chunkInd,".txt"), quote = F, row.names = F, col.names = F)
	return(NULL)
}

for (chunkInd in 1:chunkNum){
	#then combine them together
	system(paste0("cat ",output_dir,"/sample_01.umi_chunk_",chunkInd,".txt >> ",output_dir,"/sample_01.umi.txt"))
	system(paste0("rm ",output_dir,"/sample_01.umi_chunk_",chunkInd,".txt"))
}


print("Finished generating UMI barcode..")

