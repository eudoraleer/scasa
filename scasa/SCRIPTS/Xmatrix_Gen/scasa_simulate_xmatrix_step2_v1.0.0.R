#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 2
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
args = commandArgs(trailingOnly=TRUE)
whitelist <- as.character(args[1])
output_dir <- as.character(args[2])
odd_sim_fasta <- as.character(args[3])

white_list <- read.table(whitelist, header = F,stringsAsFactors=FALSE)
white_list <- white_list[,1]

txread_map=read.table(odd_sim_fasta, header = F,stringsAsFactors=FALSE)
txread_map=txread_map[,1]
txID.all=unique(txread_map)
num_reads=length(txread_map)

# UMI,12bp
umi_size <- 12
set.seed(2021)
current_umi <- do.call(paste0, replicate(umi_size, sample(c("A","T","G","C"), num_reads, TRUE), FALSE))
write.table(current_umi, paste(output_dir,"/sample_01.umi.txt",sep = ""), quote = F, row.names = F, col.names = F)
print("Finished generating UMI barcode..")

set.seed(2021)
white_list.txID.all=sample(white_list,length(txID.all),replace = FALSE)
save(txID.all,white_list.txID.all, file=paste(output_dir,"/TxCellID_mapping.RData",sep = ""))

matchID=match(txread_map,txID.all)
current_white_list=white_list.txID.all[matchID]

write.table(current_white_list, paste(output_dir,"/sample_01.white.list.txt",sep = ""), quote = F, row.names = F, col.names = F)

print("Finished generating sample white list..")
