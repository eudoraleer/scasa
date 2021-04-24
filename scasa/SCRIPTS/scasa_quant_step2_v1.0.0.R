#!/usr/bin/Rscript
##############################################################################
#                  SCASA: QUANTIFICATION
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 2
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
workdir="./"
design.matrix="Xmatrix.RData"
core = 8 #default
RsourceFn="Rsource.R"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="in") workdir=as.character(res[2])
  if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
  if (res[1]=="Rsource") RsourceFn=res[2]
}

cat("\n Create_count_matrix.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n in: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
cat("\n Rsource: ",RsourceFn)
cat("\n ----------------------------------------------------- ")

source(RsourceFn)

options(stringsAsFactors=FALSE)
setwd(workdir)
load(design.matrix)

workdir=getwd()#get full path

if(!dir.exists("Ycount")) dir.create("Ycount", recursive = T)
setwd(paste(workdir,"/Ycount",sep=""))

packages <- c("data.table","foreach","doParallel")

if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages)
}

library("data.table")
library("foreach")
library("doParallel")

load(paste0(workdir,"/Sample_eqClass.RData"))
#order by Barcode, so make it faster
setorder(sample_eqclass,Barcode) #super fast
allCells=unique(sample_eqclass$Barcode)
n=length(allCells)

cat("\n split data: ")
chunkSize=100
chunkNum=round(n/chunkSize)
chunkNum=ifelse(n > chunkSize,chunkNum, core)
chunk_eqclass=list()
brPos=round(seq(1,n+1,length.out=chunkNum+1))
for (i in 1:chunkNum){
  cat(" ",i)
  start=brPos[i]
  end=brPos[i+1]-1
  myBC=allCells[start:end]
  pick=which(sample_eqclass$Barcode%in%myBC)
  myeq=sample_eqclass[pick,]
  chunk_eqclass[[i]]=myeq
  sample_eqclass=sample_eqclass[-pick,]
}
cat("\n split data: done")

#initialization for parallel computing
registerDoParallel(cores=core)

cat("\n Run crpcount_td")
cat("\n # of chunks",chunkNum)


res=foreach(id = 1:chunkNum,.combine=c) %dopar% {
   cat("\n chunk ",id)
  mychunk=chunk_eqclass[[id]]
  allCells.chunk=unique(mychunk$Barcode)

  res_chunk=lapply(allCells.chunk,function(myBC){
    myeq=mychunk[which(mychunk$Barcode==myBC),]
    mycrp=crpcount_td(myeq)
    return(mycrp)
  })
  names(res_chunk)=allCells.chunk
  #merge them together
  Y=res_chunk[[1]]
  if (length(res_chunk) > 1){
     for (i in 2:length(res_chunk)){
       Y=Map(cbind,Y,res_chunk[[i]])
     }
  }
  #assign cell barcode
  Y=lapply(Y,function(x){  
     colnames(x)=allCells.chunk
     return(x)
   })
  id=colnames(Y[[1]])[1]
  fout=paste0(workdir,"/Ycount/ycount_chunk_",id,".RData")
  save(Y,file=fout)
  return(id)
}

chunk_eqclass=NULL #release memory

cat("\n Run crpcount_td: done")


##### estimate and export to file
flist = list.files(paste(workdir,"/Ycount/", sep = ""),pattern="ycount_chunk",recursive=FALSE,full.names = TRUE)

cat("\n estimate isoform expression ...")
DM0=CCRP
DM0 <- lapply(CRP, ccrpfun,30) #run with clim of 30
isoEst=foreach(i = 1:length(flist),.combine=list) %dopar% {
  load(flist[i])
  id=colnames(Y[[1]])[1]
  cat(" ",id)
  isoform_est=estim_chunk(Y,DM0,maxiter.bt=1000)
  fout=paste0(workdir,"/Ycount/isoEst_chunk_",id,".RData")
  save(isoform_est,file=fout)
  return(NULL)
}
res=NULL #release memory

#combine them together
flist = list.files(paste(workdir,"/Ycount/", sep = ""),pattern="isoEst_chunk_",recursive=FALSE,full.names = TRUE)

cellNum=0
txNum=0
for (i in 1:length(flist)){
  load(flist[i])
  cellNum=cellNum+nrow(isoform_est)
  txNum=ncol(isoform_est)
}
isoform_count <- matrix(0, nrow = cellNum, ncol = txNum)

curPos=0
cnames=rnames=NULL
for (i in 1:length(flist)){
  load(flist[i])
  startPos=curPos+1
  n=nrow(isoform_est)
  rnames=c(rnames,rownames(isoform_est))
  cnames=colnames(isoform_est)
  endPos=startPos+n -1
  curPos=endPos
  isoform_count[startPos:endPos,]=isoform_est
}
colnames(isoform_count)=cnames
rownames(isoform_count)=rnames
isoform_count=isoform_count[order(rownames(isoform_count)),]
isoform_count=t(isoform_count)
save(isoform_count,file=paste(workdir,"/isoform_expression.RData",sep = ""))
write.table(isoform_count, paste(workdir,"/scasa_isoform_expression.txt",sep = ""), quote = F, row.names = T, sep = "\t")
cat("\n estimate isoform expression : done!")
