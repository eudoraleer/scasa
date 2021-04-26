#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 3_2
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
packages <- c("data.table","foreach","doParallel")
if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages, repos = "https://cran.r-project.org")
}

library("data.table")
library("foreach")
library("doParallel")
core=4

args = commandArgs(trailingOnly=TRUE)
fn <- as.character(args[1])
TxCellID_mapping_fn <- as.character(args[2])
core <-as.numeric(args[3])
output_dir <- as.character(args[4])
registerDoParallel(cores=core)

dat=readLines(fn)
s=as.integer(dat[1:3]);dat=dat[-c(1:3)]
txAll=dat[1:s[1]];dat=dat[-c(1:s[1])]
cbAll=dat[1:s[2]];dat=dat[-c(1:s[2])]

cat("# of eqclasses: ",s[3],"\n")
sample_eqclass=foreach(i = 1:s[3],.combine=rbind) %dopar% {
	cat(" ",i)
	library(data.table)
	one_eqclass=data.table()
	eq=unlist(strsplit(dat[i],"\t"))
	txps_size=as.integer(eq[1]);eq=eq[-1]
	tid_set=as.integer(eq[1:txps_size]);eq=eq[-c(1:txps_size)]#transcript index, zero-based index
	tid_set=tid_set+1
	count=as.integer(eq[1]);eq=eq[-1] #eqclass count
	bgroup_size=as.integer(eq[1]);eq=eq[-1] #number of cells
	bcID=NULL
	bcCount=NULL #ycount
	pick=grep("[AGTC]",eq)
	eq2=eq[-c(pick,pick+1)]

	myid=seq(1,length(eq2),2)
	bcList=as.integer(eq2[myid]) #cellbarcode index, zero-based index
	bcList=bcList+1
	bcCount=as.integer(eq2[-myid])#count
	bcInfo=cbind(bcList,bcCount)

	one_eqclass=lapply(c(1:nrow(bcInfo)),function(j){
		x=data.table(Barcode=cbAll[bcInfo[j,1]],Transcript=txAll[tid_set],Count=bcInfo[j,2],eqClass=i)
		return(x)
	})
	one_eqclass=do.call(rbind,one_eqclass)

	return(one_eqclass)
}

save(sample_eqclass,file=paste(output_dir,"/Sample_eqClass_alevin_Xmat.RData",sep = ""))

#### Get weight ####
load(TxCellID_mapping_fn)

eqClassFinal=data.table()
eqIDSet=unique(sample_eqclass$eqClass)
cat("\n # eq ",length(eqIDSet))

eqClassFinal=foreach(i = 1:length(eqIDSet),.combine=rbind) %dopar% {
	eqID=eqIDSet[i]
	myeq=sample_eqclass[sample_eqclass$eqClass==eqID,]
	count=sum(myeq$Count)
	matchID=match(myeq$Barcode,white_list.txID.all)
	myeq$bctx=txID.all[matchID]
	txweight=myeq[, sum(Count), by = bctx]
	colnames(txweight)=c("Transcript_ID","Weight")
	mytx=sort(unique(c(unique(myeq$Transcript),txweight$Transcript_ID)))
	myeq2=data.table(Transcript_ID=mytx,eqClass=eqID,Count=count,Weight=0)
	matchID=match(txweight$Transcript_ID,myeq2$Transcript_ID)
	myeq2$Weight[matchID]=txweight$Weight
	return(myeq2)
}

write.table(eqClassFinal,paste(output_dir,"/Xmatrix_eqClass.txt", sep = ""), quote = F, row.names = F,sep="\t")
