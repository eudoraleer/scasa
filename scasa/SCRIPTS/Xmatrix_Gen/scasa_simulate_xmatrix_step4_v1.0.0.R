#!/usr/bin/Rscript
##############################################################################
#                  SCASA: SIMULATE X MATRIX
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 4
#             Author: Trung Nghia Vu, Yudi Pawitan, Lu Pan
#             Last Update: 2023-30-09
##############################################################################

### Nghia/30Sep2023: 
# - fix the bug when building CRP from CRPCOUNT

### Nghia/27Feb2023: 
# - improve codes to deal with huge X-matrices and too noisy eqclasses
# - speed up the calculation
#### build X matrix

eqClassFn="Xmatrix_eqClass.txt"
H_thres=0.025
fout='Xmatrix.RData'
RsourceFn="Rsource.R"
H_thres_min=0.0
targetsize=1000

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="in") eqClassFn=res[2]
  if (res[1]=="H") H_thres=as.double(res[2])
  if (res[1]=="s") targetsize=as.double(res[2])
  if (res[1]=="out") fout=res[2]
  if (res[1]=="Rsource") RsourceFn=res[2]
  if (res[1]=="ncore") core=res[2]
}


cat("\n buildCRP.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n in: ",eqClassFn)
cat("\n H (H threshold): ",H_thres)
cat("\n s (target to break down a huge matrix with more than s isoforms): ",targetsize)
cat("\n out: ",fout)
cat("\n Rsource: ",RsourceFn)
cat("\n ncore: ",core)
cat("\n ----------------------------------------------------- ")

source(RsourceFn)

packages <- c("foreach","doParallel","data.table")

if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages, repos = "https://cran.r-project.org")
}

library(data.table)
library(foreach)
library(doParallel)
# core = detectCores()
registerDoParallel(cores=core)

#get input from eqClass.txt in the output directory of GenTC
rawmat  = read.table(eqClassFn, header=TRUE, as.is=TRUE, sep='\t')

#remove zero eqc
pick=which(rawmat$Count !=0)
rawmat=rawmat[pick,]

#remove too lowly expressed eqc (count <= 1)
pick=which(rawmat$Count >= 2)
rawmat=rawmat[pick,]

# nghia/24Feb2023: remove noisy eqclasses
txTrueCount=tapply(rawmat$Weight,rawmat$Transcript_ID, sum)
rawmat$txTrueCount=txTrueCount[rawmat$Transcript_ID]
rawmat$prop=rawmat$Weight/rawmat$txTrueCount
rawmat$propEq=rawmat$Count/rawmat$txTrueCount
minProp=tapply(rawmat$prop,rawmat$eqClass,max)
minPropEq=tapply(rawmat$propEq,rawmat$eqClass,max)
pick=which(minProp > H_thres_min | minPropEq > H_thres_min) #excluded if reads of individual transcripts in the eqclass contributes < H_thres
rawmat=rawmat[rawmat$eqClass %in% as.integer(names(pick)),]

#reindex
myeqID=unique(rawmat$eqClass)
myeqID=sort(myeqID)
matchID=match(rawmat$eqClass,myeqID)
rawmat$eqClass=matchID
######

tx2eqc = tapply(rawmat$eqClass,rawmat$Transcript,c)  ## map: tx --> eqClass 
a = sapply(tx2eqc, length)
#table(a)
eqc2tx = tapply(rawmat$Transcript, rawmat$eqClass,c) ## map: eqClass --> tx


# neighbors of each tx
fn = function(i) {eqc=as.character(tx2eqc[[i]]);
unique(unlist(eqc2tx[eqc]))
}


#system.time(NB <- foreach(i=1:length(tx2eqc)) %dopar% fn(i) )  ## about 25sec 
#names(NB) = names(tx2eqc)

NB=lapply(tx2eqc,function(x) unique(unlist(eqc2tx[x])))
          
### build transcript clusters - new codes

#convert name of tx to index
txi=seq(length(NB))
NBi=NB
names(NBi)=txi
#map from tx to crp
t2c_map=unlist(NBi,use.names=FALSE)
t2c_mapi=match(t2c_map,names(NB))
t2c_mapi_group=rep(txi, lengths(NBi))
NBi2=tapply(t2c_mapi,t2c_mapi_group,c)
NBi=NBi2

f1=sapply(NBi,function(x) min(x))
f2=sapply(NBi,function(x) min(f1[x]))

growTimes=1
isOK=FALSE
repeat{  
  f2=sapply(NBi,function(x) min(f1[x]))  
  difNum=sum(f2!=f1)  
  growTimes=growTimes+1
  f1=f2
  if (difNum==0) isOK=TRUE
  if (isOK) break();
}

f1=names(NB)[f1]
names(f1)=names(NB)

NB2=tapply(names(f1),f1,c)
NB2=lapply(NB2,function(x) sort(x))
#create TC3
TC3=NB
TC3[names(f1)]=NB2[f1]

############ 
# nghia/26Feb2023: break down huge matrices to smaller ones by gradually increasing H_thres
getTC3<-function(matIn){
    #reindex
    myeqID=unique(matIn$eqClass)
    myeqID=sort(myeqID)
    matchID=match(matIn$eqClass,myeqID)
    matIn$eqClass=matchID
    ######
    tx2eqc = tapply(matIn$eqClass,matIn$Transcript,c)  ## map: tx --> eqClass 
    a = sapply(tx2eqc, length)
    #table(a)
    eqc2tx = tapply(matIn$Transcript, matIn$eqClass,c) ## map: eqClass --> tx

    NB=lapply(tx2eqc,function(x) unique(unlist(eqc2tx[x])))

    #convert name of tx to index
    txi=seq(length(NB))
    NBi=NB
    names(NBi)=txi
    #map from tx to crp
    t2c_map=unlist(NBi,use.names=FALSE)
    t2c_mapi=match(t2c_map,names(NB))
    t2c_mapi_group=rep(txi, lengths(NBi))
    NBi2=tapply(t2c_mapi,t2c_mapi_group,c)
    NBi=NBi2

    f1=sapply(NBi,function(x) min(x))
    f2=sapply(NBi,function(x) min(f1[x]))

    growTimes=1
    isOK=FALSE
    repeat{  
      f2=sapply(NBi,function(x) min(f1[x]))  
      difNum=sum(f2!=f1)  
      growTimes=growTimes+1
      f1=f2
      if (difNum==0) isOK=TRUE
      if (isOK) break();
    }

    f1=names(NB)[f1]
    names(f1)=names(NB)

    NB2=tapply(names(f1),f1,c)
    NB2=lapply(NB2,function(x) sort(x))
    #create TC3
    TC3=NB
    TC3[names(f1)]=NB2[f1]

    return(list(TC3=TC3,NB=NB,eqc2tx=eqc2tx,tx2eqc=tx2eqc,rawmat=matIn))
  }

len_TC3=lengths(TC3)
bigTC3=which(len_TC3>targetsize)

if (length(bigTC3) > 0){
    rmEqSet=NULL
    rmIdSet=NULL
    H_thres1=H_thres_min

    rawmat$rawEq=rawmat$eqClass
    rawmat$rawId=c(1:nrow(rawmat))
    eqSize=tapply(rawmat$Weight,rawmat$eqClass, length)
    rawmat$eqSize=eqSize[rawmat$eqClass]
    rawmat1=rawmat

    repeat{
      if (length(bigTC3) == 0 | H_thres1 >= H_thres) break()
      H_thres1=H_thres1+0.005 #increase 0.005 each time
      cat("\n",H_thres1)
      mytx=names(bigTC3)
      myeq=rawmat1$rawEq[rawmat1$Transcript_ID %in% mytx]
      #myeq=rawmat1$rawEq[rawmat1$Transcript %in% mytx]
      myeq=unique(myeq)

      rawmat2=rawmat1[rawmat1$rawEq %in% myeq,]
      minProp=tapply(rawmat2$prop,rawmat2$rawEq,max)
      minPropEq=tapply(rawmat2$propEq,rawmat2$rawEq,max)

      pick=which(minProp <= H_thres1 & minPropEq <= H_thres1) #excluded if reads of individual transcripts in the eqclass contributes < H_thres
#      pick=which(minProp <= H_thres1) #excluded if reads of individual transcripts in the eqclass contributes < H_thres

      rmEqSet=c(rmEqSet,as.integer(names(pick)))

      pick=which(minProp > H_thres1 | minPropEq > H_thres1) #excluded if reads of individual transcripts in the eqclass contributes < H_thres
#      pick=which(minProp > H_thres1)
      rawmat3=rawmat2[rawmat2$rawEq %in% as.integer(names(pick)),]

      #27Feb2023/Nghia: break "weak" links by excluding high-degree isoforms (have too many neighbours) from big eqclasses (size > 50)
      rmID=NULL
      p=which(lengths(NB) > 50)
      for (i in 1:length(p)){
        pick=which(rawmat3$Transcript == names(p)[i])
        x=rawmat3[pick,]
        pick2=which(x$propEq < H_thres1 & x$Weight==0 & x$eqSize>50)
        rmID=c(rmID,pick[pick2])

        pick3=which(x$prop < H_thres1 & x$eqSize>50)
        rmID=c(rmID,pick[pick3])
      }
      if (length(rmID)>0){
        rmIdSet=c(rmIdSet,rawmat3$rawId[rmID])
        rawmat3=rawmat3[-rmID,]
      }
      cat("\n rmID ",length(rmID))
      res=getTC3(rawmat3)
      rawmat1=res$rawmat   
      len_TC3=lengths(res$TC3)
      NB=res$NB
      bigTC3=which(len_TC3>targetsize)
      cat("\n",table(len_TC3))
    }#end of repeat

  #update new rawmat
  if (length(rmEqSet)>0){
    pick=which(!rawmat$eqClass %in% rmEqSet)
    rawmat=rawmat[pick,]
  }
  if (length(rmIdSet)>0){
    pick=which(!rawmat$rawId %in% rmIdSet)
    rawmat=rawmat[pick,]
  }
  res=getTC3(rawmat)
  rawmat=res$rawmat
  NB=res$NB
  TC3=res$TC3
  eqc2tx=res$eqc2tx
  tx2eqc=res$tx2eqc

  }#end of if

############

OTC <- sapply(TC3, paste, collapse=' ') # pasted version
names(OTC) = names(NB)
otcmap = as.list(OTC)
names(otcmap) = names(NB)

#output: list of clusters and tx->cluster map 
clust = names(table(OTC))  # clusters
OTC = strsplit(clust,split=' ')
names(OTC) = clust

#get CRP count
CRPCOUNT=list()
#system.time(
#CRPCOUNT <- foreach(i=1:length(OTC)) %dopar%{
for (i in 1:length(OTC)){
  if (i %% 100 == 0) cat(" ",i)
  #myOTC=unlist(OTC[i])
  myOTC=sort(OTC[[i]])
  names(myOTC)=NULL
  #get binary codes of eqc
  #myeqc=unique(unlist(sapply(myOTC,function(x) tx2eqc[x])))
  myeqc=unique(unlist(tx2eqc[myOTC]))
  #eqc2tx[myeqc]
  bcode=lapply(eqc2tx[myeqc],function(x) as.integer(!is.na(match(myOTC,x))))
  #bcode

  #get corresponding count - new codes
  pick=which(rawmat$eqClass %in% myeqc)

#  countname=paste(rawmat$Transcript[pick],rawmat$eqClass[pick],sep="__")
#  count=rawmat$Weight[pick] 
#  myname=expand.grid(myOTC,myeqc)
#  myname=paste(myname[,1],myname[,2],sep="__")
#  myCount=rep(0,length(myname))
#  matchID=match(countname,myname)
#  myCount[matchID]=count

  #create TC count matrix
#  mycrpCount=matrix(unlist(myCount),nrow=length(bcode),ncol=length(myOTC),byrow=TRUE)
#  colnames(mycrpCount)=myOTC
#  rownames(mycrpCount)=sapply(bcode, function(x) paste(x,collapse=""))

  
  #27Feb2023/Nghia: faster with data.table
  x=data.table(rawmat[pick,])
  y=dcast(x,eqClass ~ Transcript_ID,fill=0,value.var="Weight")
  y=y[match(names(bcode),y$eqClass),]
  #rownames(y)=sapply(bcode, function(x) paste(x,collapse=""))
  y[,1:=NULL]
  mycrpCount=as.matrix(y)
  rownames(mycrpCount)=sapply(bcode, function(x) paste(x,collapse=""))


  #sum-up if bcodes of rows are the same - 07 Aug 2020 Nghia:
  pick=which(duplicated(rownames(mycrpCount)))
  if (length(pick)>0){
    mycrpCount1=mycrpCount[-pick,drop=FALSE,]
    mycrpCount2=mycrpCount[pick,drop=FALSE,]
    for (myi in 1:nrow(mycrpCount2)){
      myj=rownames(mycrpCount1) %in% rownames(mycrpCount2)[myi]
      mycrpCount1[myj,]=mycrpCount1[myj,]+mycrpCount2[myi,]
    }
    mycrpCount=mycrpCount1
  }
  
  if (nrow(mycrpCount)>1){ # if there are more than 1 row
    mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
    #check if duplicated rownames
#    repID=table(rownames(mycrpCount))
#    repID=repID[which(repID>1)]
    repID=which(duplicated(rownames(mycrpCount)))
    names(repID)=rownames(mycrpCount)[repID]
    if (length(repID)>0){
      rmID=which(rownames(mycrpCount) %in% names(repID))
      sumcrpCount=NULL
      for (j in 1:length(repID)){
        sumcrpCount=rbind(sumcrpCount,colSums(mycrpCount[rownames(mycrpCount) %in% names(repID)[j],]))
      }
      rownames(sumcrpCount)=names(repID)
      mycrpCount=mycrpCount[-rmID,]
      mycrpCount=rbind(mycrpCount,sumcrpCount)
      mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
      mycrpCount=as.matrix(mycrpCount)
    }
  }
  
#  return(mycrpCount)
  CRPCOUNT[[i]]=mycrpCount
  }
#)

#### now compute CRP, use H_thres (default=0) to filter out too low proportion sharing between two transcripts
cat("\n Computing CRP ... ")
CRP=list()
for (k in 1:length(CRPCOUNT)){
  if (k %% 100 == 0) cat(" ",k)
  mycrpCount=CRPCOUNT[[k]]
  
  #fix the bug when colSums(mycrpCount)==0 
  txSum=colSums(mycrpCount)

  ##30Sep2023/Nghia: Fix bugs when some columns are excluded - this rarely happens
  if(length(txSum>0)==0) next();

  if (length(txSum==0)>0){
    #exclude ones with txSum==0
    p=which(txSum>0)
    #update matrix
    mycrpCount=mycrpCount[,p,drop=FALSE]
    #update bcode
    bcode=sapply(rownames(mycrpCount),function(x){
      x1=strsplit(x,"")[[1]]
      x2=x1[p]
      x3=paste(x2,collapse="")
      return(x3)
    })
    names(mycrpCount)=NULL
    rownames(mycrpCount)=bcode

    #check if duplicated rownames
    repID=which(duplicated(rownames(mycrpCount)))
    names(repID)=rownames(mycrpCount)[repID]
    if (length(repID)>0){
      rmID=which(rownames(mycrpCount) %in% names(repID))
      sumcrpCount=NULL
      for (j in 1:length(repID)){
        sumcrpCount=rbind(sumcrpCount,colSums(mycrpCount[rownames(mycrpCount) %in% names(repID)[j],]))
      }
      rownames(sumcrpCount)=names(repID)
      mycrpCount=mycrpCount[-rmID,]
      mycrpCount=rbind(mycrpCount,sumcrpCount)
      mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
      mycrpCount=as.matrix(mycrpCount)
    }

  }

  #if(ncol(mycrpCount)==0) next(); 
  
  #remove eqc with too small proportions - this is used for the clustering, too small proportions should not be consisdered as an real connection between two/more transcripts
  y=t(t(mycrpCount)/colSums(mycrpCount))
  y=ifelse(y<= H_thres,0,y) #suppress to zero
  
  z=apply(y,1,max)
  pick=z>H_thres
  x1=mycrpCount[pick,drop=FALSE,]
  
  mycrpCount_filtered=x1
  #decode eq in x1
  myclust=rep(1,ncol(x1))
  
  ### get this back
  myclust=c(1:ncol(x1))
  for (i in 1:nrow(x1)){
    x2=which(x1[i,] >0)
    x3=which(myclust %in% myclust[x2])
    myclust[x3]=min(myclust[x3])
  }
  
  #generate new CRP
  clustID=as.integer(names(table(myclust)))
  #newx=list()
  for (i in 1:length(clustID)){  
    pick=which(myclust == clustID[i])
    X=mycrpCount_filtered[,pick,drop=FALSE]
    pick.row=rowSums(X)>0
    X=X[pick.row,drop=FALSE,]
    
    #collect the bcodes
    bcode=lapply(rownames(X), function(x) as.integer(unlist(strsplit(x,""))))
    bcode=do.call(rbind,bcode)    
    bcode.max=apply(bcode,2,max)
    pick.col=which(bcode.max>0)
    
    bcode.X=bcode[,pick.col,drop=FALSE]
    bcode.X.names=apply(bcode.X,1,paste0,collapse="")
    
    X_names=colnames(mycrpCount_filtered)[pick.col]
    X_names=paste(X_names,collapse=" ")
    
    X=mycrpCount_filtered[pick.row,pick.col,drop=FALSE]
    rownames(X)=bcode.X.names
    
    #    #normalise to get crp
    X=t(t(X)/colSums(X))
    X=ifelse(is.na(X),0,X)
    #add up results
    newx=list()
    newx[[X_names]]=X
    CRP=c(CRP,newx)
  }
  
}

CCRP <- lapply(CRP, ccrpfun,30) #run with clim of 30

save(CRP,CCRP,CRPCOUNT,file=fout)
