#!/usr/bin/Rscript
##############################################################################
#                  SCASA: GENERATE EQCLASS
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 1_1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################

args = commandArgs(trailingOnly=TRUE)
fn <- as.character(args[1])
output_dir <- as.character(args[2])
core <- as.numeric(args[3])

dat=readLines(fn)

s=as.integer(dat[1:3]);dat=dat[-c(1:3)]
txAll=dat[1:s[1]];dat=dat[-c(1:s[1])]
cbAll=dat[1:s[2]];dat=dat[-c(1:s[2])]

#now read individual eqclasses
cat("# of eqclasses: ",s[3],"\n")

packages <- c("foreach","doParallel","data.table")

if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages)
}

library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=core)

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
	bcCount=NULL #ycount of model
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

sample_eqclass$Weight=1

save(sample_eqclass,file=paste(output_dir,"/Sample_eqClass.RData", sep = ""))
