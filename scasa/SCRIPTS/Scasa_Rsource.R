#28Mar2021/Nghia: revise the estim_chunk() function to deal with CRP of two isoforms with a poor condition
#23Mar2021/Nghia: function crpcount_td(): need an extra step to check if these eqclasses should be matched here or not, otherwise, it can be counted more than 1 if one tx can belong to 2 eqclasses 
### 27Dec2020/Nghia: function crpcount_td(): add the weight of tx in eqclass
## 16 Nov 2020 / Lu:
# Revised crpcount()
# Added new function hamming_correction()
## 25 Jul 2020 / Nghia:
# revise function ccrpfun() to obtain stable CCRPs by improve the convergence of k-mean 
## 11 Jul 2020 / Nghia:
# revised old and slow crpcount function for 10x genomics data, will improve speed later
## 27 Nov 2019 / Nghia:
# add function getSE() to compute the standard error for the estimates
## Nghia 05 May 2019: 
# - crpcount: vectorisation for a extremely fast computation
## YP 23 April 2018!!
##    revise tcfun(.) not separating the NB1 and removing from the rest,
##    which generated wrong singletons    
##
## YP 25 April 2016: add EMM = modified EM
## 13 Dec 2016: adding estimation by LSP
## 7 Dec 2016: marking the paralogs
## 5 Dec 2016: parallel version of the 4Dec version

## Processing transcript-response profile (TRP), to make
## transcription cluster (TC) and cluster response profile (CRP)
##
## tcfun function:
## input TRP  = the list of trp matrices (must be matrix!!)
## output = transcription clusters
##
## crpfun function:
## input TRP 
##   computes transcription clusters using tcfun
## output= cluster-response profile


## ................................................
## ................................................ tcfun: cluster growing function
## ................................................
tcfun = function(TRP, lim=0.025)
{
# setup neighborhood (NB) list from TRP
# drop likely false neighbors from the TRP if tx-count leakage is too low
NB.n= lapply(TRP,apply,2,sum)  # neighbors with strength=n count
fn = function(x, lim=lim) {  ## keep neighbors with leakage rate >lim
     x = x/max(x); 
     return(unique(names(x[x>lim])))
  }
NB = lapply(NB.n,fn,lim=lim)

## checking singletons, then remove them from the rest
#NBvec= sapply(NB,paste, collapse=' ')
#tst1 = names(NB) == NBvec
#  NB1 = NB[tst1]
#NB2.tmp = NB[!tst1] 
#  fn = function(i) setdiff(NB2.tmp[[i]],NB1)
#  NB2 = foreach(i= 1:length(NB2.tmp)) %dopar% fn(i)  # parallel
#    names(NB2) = names(NB2.tmp)

# TC = union of neighbors of TRPs 
## TC2 = cluster growing from NB
growfn = function(i,NB) sort(unique(unlist(NB[NB[[i]]])))  ## grow from neighborhood NB
system.time(TC2 <- foreach(i= 1:length(NB)) %dopar% growfn(i, NB))
  names(TC2) = names(NB)

# TC3 = grow from TC2
system.time(TC3 <- foreach(i= 1:length(NB)) %dopar% growfn(i, TC2))
TC =  sapply(TC3, paste, collapse=' ') # pasted version, needed for tabulating
  names(TC) = names(NB)
  tcmap = TC  ## map of tx --> TC

#output: list of clusters and tx->cluster map 
clust = names(table(TC))  # clusters
TC = strsplit(clust,split=' ')
  names(TC) = clust

return(list(TC=TC, tcmap = tcmap))
}

## .......................................................
## .......................................................crpfun:  Cluster-response profile
## .......................................................
## input: TRP and TC=tcfun(TRP)$TC
## output: CRP
crpfun = function(TRP, TC)
{
## ALL transcripts
 TX = names(TRP)  
## TC = list of clusters
 clust.size = sapply(TC,length)  ##

### singleton clusters
 C1 = which(clust.size==1)   
 fn = function(x){a= matrix(1); colnames(a) = x ; rownames(a)='1'; return(a)}
 TC1 = lapply(names(C1),fn)
   names(TC1) = names(C1)

## clusters with >1 transcripts:
  TC2 = TC[clust.size>1]
   nc = length(TC2)

## .............upat =  unique occupancy pattern 
## assigning unique patterns based on a fixed-set of transcripts TXc typically from a cluster
## ADD ZERO cols to transcripts missing in trp mat. 
##
upat = function(trp, TXc){
  TX0 = setdiff(TXc, colnames(trp))  # missing tx in the trp
  n0= length(TX0)
  if (n0>0){
     trp0= matrix(0,nrow=nrow(trp),ncol=n0);
     colnames(trp0) = TX0
     trp = cbind(trp, trp0)
  }
  out =  trp[,TXc, drop=FALSE]  ## keep matrix property even for vectors
  out01= ifelse(out>0,1,0)      ## 0-1 pattern: 0=no expression. NOTE: Think about FALSE positive!!
  pat = apply(out01,1,paste,collapse='')
  rownames(out) = pat

  ## collapsing row to unique lines: be careful with 1 unique pattern!!
  nu = length(unique(pat))   ## number of unique pattern
  if (nu>1) out = apply(out,2,tapply,rownames(out),sum)
  if (nu==1) {
      a = matrix(apply(out,2,sum),ncol=ncol(out))
      colnames(a)= colnames(out)
      rownames(a) = unique(pat)
      out = a
  }
  return(out)
}

## for parallel
 incrpfun = function(i){
   TXc = TC2[[i]]  ## set of transcripts in cluster i
     tmp = TXc %in% TX   ## if TX is not complete, TXc may contain external tx
     TXc = TXc[tmp] 
     nt = length(TXc)
   TRPc= TRP[TXc]
   ordtrp = lapply(TRPc, upat, TXc)  ## apply upat to each trp 
   ## all unique patterns
   allpat = unique(unlist(lapply(ordtrp, rownames)))

   ## setup crp for a cluster
   crp = matrix(0, nrow=length(allpat), ncol=nt)
     rownames(crp) = allpat
     colnames(crp) = TXc
   for (nam in TXc){
     rname = rownames(ordtrp[[nam]])
     crp[rname,nam] = ordtrp[[nam]][,nam]
   }
  ## normalize to sum to 1:
  crp = t(t(crp)/colSums(crp))
  return(crp)
}
CRP = foreach(i= 1:length(TC2)) %dopar% incrpfun(i)
names(CRP) = names(TC2)

return(c(TC1,CRP)) ## combine the singletons and others 
}

## ...............................................
## ............................................... ccrpfun: combining paralogs (similar transcripts)
## ...............................................          within each crp
## usage: 
## set.seed(2016)  ## to fix the kmeans clustering results 
## CCRP <- lapply(CRP, ccrpfun, clim=50)
## clim= minimum condition number
##
ccrpfun = function(x, clim=50) {
 set.seed(2016)  ## to fix the kmeans clustering results
 #remove the transcript without contribution to the CRP
 pick=which(colSums(x)>0)
 x=x[,pick,drop=FALSE]
 newx = x
 repeat{
   sval = svd(newx)$d;     ## singular values
   #con = max(sval)/sval    ## condition numbers
   con = abs(max(sval)/sval)    ## condition numbers; note from Nghia: might use abs() to avoid the case of -Inf when dividing to zero. This happened randomly, hard to reproduce
   nsc = sum(con< clim)    ## number of clearly separated transcripts
   ## break if OK
   if (nsc == ncol(newx)) break 
   #27Dec2020/Nghia: use <= instead of ==, to avoid the colSums values less than 1

   ## but collapse if some tx are too similar, using kmeans clustering
   ##   to combine the similar tx's
   ## 24Jul2020/Nghia: increase the iterations and the number of randome sets to get converged results 
   subc = kmeans(t(x),nsc,iter.max = 100000,nstart=100) ## always start with original tx from the crp! 
   newx = t(subc$center)
   clust = subc$cluster
     clust.name = tapply(names(clust), clust, paste, collapse=' ')
     colnames(newx)= clust.name
 }
 return(newx)
}

## ....................................  LS-Pos function
## ....................................  estimation using positive-constrained LS
##  
#require(limSolve)
lsp = function(X,y){
 A = X
 B = y[rownames(X)] ## equalize order of y = rows of X
 p = ncol(X)
 G <- diag(p)
 H <- rep(0,p)
 lsp = lsei(A = A, B = B, G = G, H = H, type=2)
 bt= lsp$X
 err = y- c(X %*% bt)
 L2 = sum(err^2)
 return(list(beta=bt, error=err, L2 = L2))
}

##  
## ............. estimate using EM algorithm
## ............. beta0 = starting value 
##
EM = function(X,y,beta0=rep(sum(y)/ncol(X), ncol(X)), maxiter=50, maxerr=0.01, lim=0.01){
  csum = colSums(X)  ## should be =1, but here just in case
  for (i in 1:maxiter){
    X1 = t(t(X)*beta0)                 ## expected per cell under current est
      rsum = rowSums(X1)
      X2 = X1/ifelse(rsum>0, rsum,1)   ## normalize each row to get probabilities; avoid div by zero
    beta = c(t(X2) %*% y) /csum        ## update:add up info across rows, then normalize by colSum
    err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim)
    if (max(err) < maxerr) break
    beta0 = beta
  }
  return(beta)
}

##
## modify=TRUE: modify beta with digamma function
##
EMM = function(X,y,beta0=rep(sum(y)/ncol(X), ncol(X)), maxiter=100, 
               maxerr=0.01, lim=0.01, modify=TRUE){
   csum = colSums(X)  ## should be =1, but here just in case
                      ## sailfish X version, is row-normalized, so does not add-up to one
   logNorm = digamma(sum(y))
   for (i in 1:maxiter){
    gamma0= beta0
    if (modify) gamma0=exp(digamma(beta0+0.00000001) - logNorm)
    # beta0=ifelse(is.na(beta0),0,beta0)
    #print(c(i,beta0, sum(beta0)))
    X1 = t(t(X)*gamma0)                 ## expected per cell under current est
    rsum = rowSums(X1)
    X2 = X1/ifelse(rsum>0, rsum,1)    ## normalize each row to get probabilities; avoid div by zero
    beta = c(t(X2) %*% y) # /csum        ## update:add up info across rows, then normalize by colSum
    err = abs(beta-beta0)/ifelse(beta0>lim, beta0, lim)
    if (max(err) < maxerr) break
    beta0 = beta
    }
   return(beta)
 }

###crpcount function - for fast computation -but not suitable for 10x genomics
#crpcount = function(filename, weight_thres=0.9)
#{
#  options(stringsAsFactors=FALSE)
#  rawcount = read.table(filename,header=T)
#  sample.name1 = strsplit(filename,"/")[[1]]
#  sample.name = sample.name1[length(sample.name1)-1]
#
#  ## defining eqClasses in the CRPs
#  ## all CRPs are assumed matrices -- even 1x1 singletons
#  ##  *** Tx columns within CRP are sorted***
#  crp.eqfun = function(crp){
#    tx = colnames(crp)
#    txmat = outer(rep('',nrow(crp)), tx, paste0)
#    upat = matrix('', nrow=nrow(crp), ncol=ncol(crp))
#    upat[crp>0] = txmat[crp>0]
#    out = apply(upat,1,paste, collapse='')
#    return(out)
#  }
#  ## full list eqClasses from CRPs: 
#  CRP.eq <- lapply(CRP, crp.eqfun)# 15 Otc 2019/Nghia: use lapply instead of sapply, this avoids error when having only 1 CRP
#  CRP.names = rep(names(CRP.eq), lengths(CRP.eq))
#  CRP.eqID = paste(CRP.names, unlist(CRP.eq)) 
#  ## check: each eqClass is unique, no repeat
#  ## length(unique(unlist(CRP.eq))); length(unlist(CRP.eq))
  #
#  #map from crp to tx
#  c2t=sapply(names(CRP), strsplit, split=' ')
#  cLen=lengths(c2t)
#  cVec=rep(names(cLen),cLen)
#  tVec=unlist(c2t)  
#  names(cVec) = tVec
#  rawcount$crp = cVec[rawcount$Transcript]
  #
#  #do the weight filter here
#  matchID=match(rawcount$crp,names(cLen))
#  rawcount$eqlen=cLen[matchID]  
#  pick = (rawcount$eqlen>1)| 
#    (rawcount$eqlen==1 & rawcount$Weight>weight_thres )
#  # pick = rep(TRUE, nrow(rawcount))
  #
#  rawcount$eqCrp=paste(rawcount$crp,rawcount$eqClass)
#  eqClass = tapply(rawcount$Transcript[pick],rawcount$eqCrp[pick],
#                   function(x){paste(sort(x),collapse='')})
#  eqCount=tapply(rawcount$Count[pick],rawcount$eqCrp[pick],min)
#  eqCRP=tapply(rawcount$crp[pick],rawcount$eqCrp[pick],unique)
#  eqID = paste(eqCRP, eqClass)
#  eqCount2=tapply(eqCount,eqID,sum)
#  eqID2=names(eqCount2)
  #
#  ##  .. transfer data to CRP objects
#  ycount = rep(0, length(CRP.eqID))
#  names(ycount) = CRP.eqID
#  pick = eqID2 %in% CRP.eqID
#  ycount[eqID2[pick]]= eqCount2[pick]
  #
#  ## repack into CRPs
#  CRP.y = tapply(ycount, CRP.names, c)
#  out = Map(cbind, CRP, sample1=CRP.y[names(CRP)])
  #
#return(list(crpcount=out,samplename=sample.name))
#}


crpcount = function(filename)
{
 options(stringsAsFactors=FALSE)
 rawcount = read.table(filename,header=T)
 #sample.name1 = strsplit(filename,"/")[[1]]
 #sample.name = sample.name1[length(sample.name1)-1]
 #sample.name=strsplit(sample.name,"_")[[1]][split2]
 if(length(grep(".*\\/eqClass.txt",filename))>0){
filename = gsub("\\/\\/","\\/",filename)
if(length(as.character(unlist(strsplit(filename,split="\\/")))) > 2){
sample.name = gsub(".*\\/(.*)\\/eqClass.txt","\\1",filename)
}else{
sample.name = gsub("(.*)\\/eqClass.txt","\\1",filename)
}
}else{
sample.name = gsub(".*\\/(.*)_eqClass.txt","\\1",filename)
sample.name = gsub("(.*)_eqClass.txt","\\1",sample.name)

} ###
 #remove zero eqc
 pick=which(rawcount$Count !=0)
 rawmat=rawcount[pick,]
 #reindex
 myeqID=unique(rawcount$eqClass)
 myeqID=sort(myeqID)
 matchID=match(rawcount$eqClass,myeqID)
 rawcount$eqClass=matchID
 ### 27Dec2020/Nghia
 #set weight for eqclass
 eqWeight=table(rawcount$eqClass)
 eqWeight=1/eqWeight
 rawcount$Weight=eqWeight[match(rawcount$eqClass,as.integer(names(eqWeight)))]
 ###
 tx1 = unique(rawcount$Transcript)
 nID=c(1:length(CRP))
 rapmap_crp=NULL
 extra_pattern = NULL
 for (j in nID)
 {
 #cat(j,'\n')
 crp = CRP[[j]]
 if(length(crp)==1)
 {
 tx = colnames(crp)
 sample1 = sum(rawcount[rawcount$Transcript==tx & rawcount$Weight>0.9,"Count"])
 crp1 = cbind(crp,sample1)
 crp.y = crp1
 }else{
   tx = colnames(crp)
   #27Dec2020/Nghia
   pick=colSums(crp)>0
   tx1=tx[pick] #consider only eqc that contribute to the CRP
   #
   raw.c =NULL
   pat=NULL
   for(i in 1:length(tx1))
      raw.c = rbind(raw.c,rawcount[rawcount$Transcript==tx1[i],])
   eq = unique(raw.c$eqClass)
   if(length(eq)>0)
 {
 raw2 = matrix(0,length(eq),length(tx))
 colnames(raw2) = tx
 rownames(raw2) = eq
 for(z in 1:nrow(raw.c))
 raw2[as.character(raw.c[z,grep("eqClass", colnames(raw.c), ignore.case = TRUE )]),raw.c[z,1]] = raw.c[z,grep("Count", colnames(raw.c), ignore.case = TRUE )]##3 is counts, 6 is eqclass

 raw3 = ifelse(raw2>0,1,0)
 pat = apply(raw3,1,function(x)paste(x,collapse=""))
 raw3 = data.frame(pat=pat,raw3)
 raw2.1 = data.frame(pat=pat,raw2)
 a1 = apply(raw2,1,max)
 a2 = tapply(a1,pat,sum)
 pat=unique(pat)
 raw4 = data.frame(pat=names(a2),sample1=a2)
 crp1 = data.frame(pat=rownames(crp),crp)

 raw4 <- hamming_correction(raw4, crp1, hamming_dist = 1)
 
 m = merge(crp1,raw4,by="pat",all.x=TRUE,sort=FALSE)
 rownames(m) = as.character(m$pat)
 m[is.na(m)] = 0
 m = m[,grep("pat", colnames(m), ignore.case = TRUE, invert = TRUE)]
 crp.y = m[rownames(crp),]
 
 }else{
 crp.y = data.frame(crp,sample1=0)
}

}
crp.y=as.matrix(crp.y)
if(j ==1) rapmap_crp = list(crp.y)
if(j >1) rapmap_crp[[j]] = crp.y
}
names(rapmap_crp) = names(CRP)
return(list(crpcount=rapmap_crp,samplename=sample.name))
}

#crpcount_td = function(Barcode,Transcript,Count,eqClass,Weight)
crpcount_td = function(mytd)
{
 options(stringsAsFactors=FALSE) #must have
# rawcount=data.table(Transcript=Transcript,Count=Count,eqClass=eqClass,Weight=Weight)
# sample.name=Barcode

 #map from crp to tx
c2t=sapply(names(CRP), strsplit, split=' ') # pasted version
#map from tx to crp
t2c=unlist(c2t,use.names=FALSE)
names(t2c)=rep(names(c2t), lengths(c2t))
t2c=tapply(names(t2c),t2c,c)


 rawcount=data.frame(mytd[,-c("Barcode")],stringsAsFactors=FALSE)
 rawcount$Transcript=as.character(rawcount$Transcript)
 sample.name=mytd$Barcode[1]

 #remove zero eqc
 pick=which(rawcount$Count !=0)
 rawmat=rawcount[pick,]
 #reindex
 myeqID=unique(rawcount$eqClass)
 myeqID=sort(myeqID)
 matchID=match(rawcount$eqClass,myeqID)
 rawcount$eqClass=matchID
 ### 27Dec2020/Nghia
 #set weight for eqclass
 eqWeight=table(rawcount$eqClass)
 eqWeight=1/eqWeight
 rawcount$Weight=eqWeight[match(rawcount$eqClass,as.integer(names(eqWeight)))]

 ###
 tx1 = unique(rawcount$Transcript)
 nID=c(1:length(CRP))
 rapmap_crp=NULL
 extra_pattern = NULL
 for (j in nID)
 {
 #cat(j,'\n')
 crp = CRP[[j]]
 if(length(crp)==1)
 {
 tx = colnames(crp)
 sample1=0
 pick=which(rawcount$Transcript %in% tx & rawcount$Weight>0.9)
 if (length(pick)>0) sample1 = sum(rawcount[pick,"Count"])
 crp1 = cbind(crp,sample1)
 crp.y = crp1
 }else{
   tx = colnames(crp)
   #27Dec2020/Nghia
   pick=colSums(crp)>0
   tx1=tx[pick] #consider only eqc that contribute to the CRP
   #
   raw.c =NULL
   pat=NULL
   for(i in 1:length(tx1))
      raw.c = rbind(raw.c,rawcount[rawcount$Transcript==tx1[i],])
   eq = unique(raw.c$eqClass)
   #23Mar2021/Nghia: need an extra step to check if these eqclasses should be matched here or not, otherwise, it can be counted more than 1 if one tx can belong to 2 eqclasses   
if(length(eq)>0){
  currentCRPname=names(CRP[j])
  myRawcount=rawcount[rawcount$eqClass %in% eq,]
  rmList=NULL
  for (i in 1:length(eq)){
    mytx=myRawcount$Transcript[myRawcount$eqClass==eq[i]]
    mytx=sort(mytx)
    eq_candidates=sort(unique(unlist(t2c[mytx])))
    
    #keep only relevant eq
    keepID=NULL
    for (k in 1:length(eq_candidates)){
      acrp=CRP[[eq_candidates[k]]]
      acrp_tx=colnames(acrp)[colSums(acrp)>0]
      if (sum(acrp_tx %in% mytx)>0) keepID=c(keepID,k)
    }
    eq_candidates=eq_candidates[keepID]

    #compare the set of tx of the eqclass and the crp
    myscore=NULL
    for (k in 1:length(eq_candidates)){
      #27Mar2021:Nghia: tricky problem, some tx in CRP are not expressed, so should not be included
      x=c2t[[eq_candidates[k]]]
#      acrp=CRP[[eq_candidates[k]]]
#      x=colnames(acrp)[colSums(acrp)>0]
      sc=sum(x %in% mytx)/length(unique(c(x,mytx)))
      myscore=c(myscore,sc)
    }
    bestEq=eq_candidates[which.max(myscore)]
    if (currentCRPname!=bestEq) rmList=c(rmList,i)
  }
  if (length(rmList)>0) eq=eq[-rmList]
  if (length(eq)>0){ #update raw.c
  raw.c=raw.c[raw.c$eqClass %in% eq,drop=FALSE,]
  }
}

if(length(eq)>0)
 {
 raw2 = matrix(0,length(eq),length(tx))
 colnames(raw2) = tx
 rownames(raw2) = eq
 for(z in 1:nrow(raw.c)){
 d1=as.character(raw.c[z,grep("eqClass", colnames(raw.c), ignore.case = TRUE )])
 d2=raw.c[z,1]
 raw2[d1,d2] = raw.c[z,grep("Count", colnames(raw.c), ignore.case = TRUE )]
}
 raw3 = ifelse(raw2>0,1,0)
 pat = apply(raw3,1,function(x)paste(x,collapse=""))
 raw3 = data.frame(pat=pat,raw3)
 raw2.1 = data.frame(pat=pat,raw2)
 a1 = apply(raw2,1,max)
 a2 = tapply(a1,pat,sum)
 pat=unique(pat)
 raw4 = data.frame(pat=names(a2),sample1=a2)
 crp1 = data.frame(pat=rownames(crp),crp)

 raw4 <- hamming_correction(raw4, crp1, hamming_dist = 1)
 
 m = merge(crp1,raw4,by="pat",all.x=TRUE,sort=FALSE)
 rownames(m) = as.character(m$pat)
 m[is.na(m)] = 0
 m = m[,grep("pat", colnames(m), ignore.case = TRUE, invert = TRUE)]
 crp.y = m[rownames(crp),]
 
 }else{
 crp.y = data.frame(crp,sample1=0)
}

}
crp.y=as.matrix(crp.y)
pick=which(colnames(crp.y)=="sample1")
#colnames(crp.y)[pick]=sample.name
crp.y=crp.y[,pick]
names(crp.y)=NULL
if(j ==1) rapmap_crp = list(crp.y)
if(j >1) rapmap_crp[[j]] = crp.y
}
names(rapmap_crp) = names(CRP)
return(rapmap_crp)
}


# 2020-11-16:
# Function added by: Lu, to fix the problems with no comparable equivalence class
hamming_correction <- function(raw_data, crp_data, hamming_dist = 1){   # by default, hamming correction of distance = 1
  raw_current=lapply(raw_data$pat, function(x) as.integer(unlist(strsplit(x,""))))
  crp_current=lapply(crp_data$pat, function(x) as.integer(unlist(strsplit(x,""))))
  
  #compute hamming distance
  res=lapply(raw_current, function(x){
    z=sapply(crp_current,function(y) sum(x!=y))
    return(z) 
  })
  res=do.call(rbind,res)
  dim(res)
  rownames(res)=raw_data$pat
  colnames(res)=crp_data$pat
  res=as.matrix(res)
  compare_hamming=res

#  raw_current <- data.frame(t(data.frame(strsplit(raw_data$pat, split = ""))))
#  for(i in 1:ncol(raw_current)){
#    raw_current[,i] <- as.numeric(as.character(raw_current[,i]))
#  }
#  row.names(raw_current) <- raw_data$pat
  #
#  crp_current <- data.frame(t(data.frame(strsplit(crp_data$pat, split = ""))))
#  for(i in 1:ncol(crp_current)){
#    crp_current[,i] <- as.numeric(as.character(crp_current[,i]))
#  }
#  row.names(crp_current) <- crp_data$pat
  #
#  compare_hamming <- NULL
#  for(p in 1:nrow(raw_current)){
#    temp <- NULL
#    for(q in 1:nrow(crp_current)){
#      current <- length(which(abs(crp_current[q,] - raw_current[p,]) > 0))
#      temp <- cbind(temp, current)
#    }
#    compare_hamming <- rbind(compare_hamming, temp)
#  }
  #
#  compare_hamming <- data.frame(compare_hamming)
#  colnames(compare_hamming) <- row.names(crp_current)
#  row.names(compare_hamming) <- row.names(raw_current)

  raw_corrected <- raw_data
  # current <- apply(compare_hamming,1,function(x){which(x == 1)})
  for(i in 1:nrow(compare_hamming)){

    #if this eqclass is found in crp
    if(sum(compare_hamming[i,] == 0) > 0){ 
      raw_corrected$pat[i] <- colnames(compare_hamming)[which(compare_hamming[i,] == 0)]
    }else{
      closeNum=sum(compare_hamming[i,] <= hamming_dist) #number of eqclasses with hamming_dist less than the threshold
      if( closeNum > 1){ #more than one close-eqclass
        # Lu noted here: If we choose only the one with highest value, then might be a bias towards high value in real data
        # might affect those with real cases of low values
        current <- crp_data[which(compare_hamming[i,] <= hamming_dist),]
        raw_corrected$pat[i] <- row.names(current)[which.max(apply(current,1,function(x){median(as.numeric(as.character(x[x>0])))}))]
      }

      if(closeNum == 1){ #only one-close eqclass
        raw_corrected$pat[i] <- colnames(compare_hamming)[which(compare_hamming[i,] <= hamming_dist)]
      }

      if(closeNum == 0){ #no close-eqclass
        current <- crp_data[which(compare_hamming[i,] <= min(compare_hamming[i,])),]
        raw_corrected$pat[i] <- row.names(current)[which.max(apply(current,1,function(x){median(as.numeric(as.character(x[x>0])))}))]
      }
    }
  }
  
#  if(sum(duplicated(raw_corrected$pat)) > 0){
  #
  #current <- raw_corrected[!(!(duplicated(raw_corrected$pat) | duplicated(raw_corrected$pat, fromLast = TRUE))),]
  #raw_corrected <- raw_corrected[(!(duplicated(raw_corrected$pat) | duplicated(raw_corrected$pat, fromLast = TRUE))),]
#  current <- split(current, current$pat)
  ## 2020-11-17: Modified by Lu, changed from select from the one that has max median to include all
  ## current <- lapply(current, function(x){x <- x[which.max(apply(data.frame(x[,grep("pat", colnames(x), ignore.case = TRUE, invert = TRUE)]),1,median)),]})
  #current <- lapply(current, function(x){x <- data.frame(pat = unique(x$pat), sample1 = sum(x$sample1))})
#  current <- do.call(rbind.data.frame, current)
#  raw_corrected <- data.frame(rbind(raw_corrected, current))
#  }
  
   a = tapply(raw_corrected$sample1,raw_corrected$pat,sum)
   raw_corrected = data.frame(pat=names(a),sample1=a)
 
  row.names(raw_corrected) <- raw_corrected$pat
  return(raw_corrected)
}

## 25 April 2017: 
## .. YP: estim.BETA and AEM allow modify=TRUE for digamma modification
##
## 18 April 2017
## ..  estim.BETA does not use big matrix Z  
##
## joint estimation of A and theta for sequgio model
## SEE sequgio-model and estimation.r
##

## to use LSP= positive-constrained least-squares from library(limSolve)
#require(limSolve)

##  ......................  Alternating EM algorithm applied to one CRP
##
##   X0   = starting value for design matrix X = CRP
##   Ymat = data matrix, ncol(Ymat)= sample size n; nrow(Ymat) = num of eqClasses
##
AEM = function(X0, Ymat, maxiter.X=10, maxerr=0.01, lim=0.01,
                         maxiter.bt=100, modify=TRUE){
 X.est0 = X0  # starting design matrix
 Ysum= colSums(Ymat)
 beta0 = rep(Ysum/ncol(X0), rep(ncol(X0), length(Ysum))) # starting espression
 ## alternating estimation
 for (i in 1:maxiter.X){
   BT = estim.BETA(X.est0, Ymat, beta0=beta0, maxiter=maxiter.bt, modify=modify); 
   X.est = estim.X.em(X.est0,BT,Ymat); ## print(rbind(BT,X.est))
   ## check convergence
   err = abs(X.est-X.est0)/ifelse(abs(X.est0)>lim, X.est0, lim)  ## X.est0 could be zero
   if (max(err) < maxerr) break
   ## update for next iterate
   X.est0 = X.est
   beta0=c(t(BT))   ## NOTE: BT is a matrix, saved byrow=T
 }
return(list(X = X.est, BETA=BT))
}


##
## EM without for-loop over sample size n, only over p
##
estim.X.em = function(X.em,BETA,Ymat){
  X0 = X.em  ## starting value
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X.em)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)     ## sample size
  p = ncol(BETA)  ## number of transcripts
 BTsum = rep(colSums(BETA), rep(q,p)) ## NOTE: BTsum might be zero!!
 x1 = NULL
 for (j in 1:p) x1=cbind(x1,kronecker(BETA[,j],X.em[,j]))
 sumx1 = rowSums(x1) 
   sumx1.pos= ifelse(sumx1>0, sumx1, 0.1)
 x2 = x1/sumx1.pos
 Ycell= x2*Y
 tmp = NULL 
 for (j in 1:p){
   tmp = c(tmp,colSums(matrix(Ycell[,j], nrow=n,byrow=TRUE)))
  }
 ## NOTE BTsum might be zero!!
 ## if csum=0, return to the original X0
 X.em = matrix(tmp/ifelse(BTsum>0, BTsum, 0.1), ncol=p)
   csum = colSums(X.em)
   X.em = t(t(X.em)/ifelse(csum>0,csum,0.1))  ## normalize to sum to 1
 for (j in 1:length(csum)){if (csum[j]< 0.00001) X.em[,j]=X0[,j]}
 return(X.em)
}



##
## .. modify=TRUE: EM modified by the digamma function
## .. only EM and WITHOUT big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. Ymat size = qxn matrix: sample = columns
## .. beta0 vector of starting value for BETA: length = nxp
##
estim.BETA = function(X,Ymat,beta0, maxiter=100, maxerr=0.01, 
                         lim=0.01, modify=TRUE){  
  Y.list =  split(t(Ymat), seq(ncol(Ymat))) 
  q = nrow(X)    ## number of quasi-exons
  n = ncol(Ymat) ## number of samples
  csum = colSums(X)
  totCount= colSums(Ymat)
  logNorm = digamma(totCount+0.00000001)

  ## computing the old X2, allow modify=TRUE
  fn1 = function(bt, logNorm.i){
    gamma0= bt
    if (modify) gamma0=exp(digamma(bt+0.00000001) - logNorm.i)
    tmp = t(t(X)*gamma0)    ## NOTE: X is fixed, taken from the input
    rsum = rowSums(tmp)
    return(tmp/ifelse(rsum>0, rsum,1))
  }
  ## computing each (Xi'yi)
  fn2 = function(Xi,yi){
    bt = c(t(Xi) %*% yi) /csum  ## normalized by colSums(X)
    return(bt)
  }

for (i in 1:maxiter){
  BETA = matrix(beta0, nrow=n, byrow=TRUE)
  BT.list = split(BETA, seq(nrow(BETA))) ## list version
  X2.list = mapply(fn1, BT.list, logNorm, SIMPLIFY=FALSE)
  beta= mapply(fn2, X2.list, Y.list, SIMPLIFY=TRUE)
  beta=ifelse(is.na(beta),0,beta) #27Dec2020/Nghia
  err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
  if (max(err) < maxerr) break
    beta0 = beta
}  ## endfor i

return(matrix(beta, nrow=n, byrow=TRUE))
}


## .. older version: EM with no modification 
## .. only EM and WITHOUT big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. Ymat size = qxn matrix: sample = columns
## .. beta0 vector of starting value for BETA: length = nxp
##
estim.BETA.M0 = function(X,Ymat,beta0, maxiter=100, maxerr=0.01, lim=0.01){  
  Y.list =  split(t(Ymat), seq(ncol(Ymat))) 
  q = nrow(X)    ## number of quasi-exons
  n = ncol(Ymat) ## number of samples
  csum = colSums(X)

  ## computing the old X2 
  fn1 = function(bt){
    tmp = t(t(X)*bt)    ## NOTE: X is fixed, taken from the input
    rsum = rowSums(tmp)
    return(tmp/ifelse(rsum>0, rsum,1))
  }
  ## computing each (Xi'yi)
  fn2 = function(Xi,yi){
    bt = c(t(Xi) %*% yi) /csum  ## normalized by colSums(X)
    return(bt)
  }

for (i in 1:maxiter){
  BETA = matrix(beta0, nrow=n, byrow=TRUE)
  BT.list = split(BETA, seq(nrow(BETA))) ## list version
  X2.list = Map(fn1,BT.list)
  beta= mapply(fn2,X2.list,Y.list)
  err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
  if (max(err) < maxerr) break
    beta0 = beta
}  ## endfor i

return(matrix(beta, nrow=n, byrow=TRUE))
}


##
## .. using potentially big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. beta0 = starting value for BETA
##
estim.BETA.Z = function(X,Ymat,beta0, method='EM', maxiter=100, maxerr=0.01, lim=0.01){  
  Y = c(Ymat)
  q = nrow(X)   ## number of eqClasses
  n = ncol(Ymat)
  Z = diag(n) %x% X   ## kronecker product
 if (method=='OLS') {
    ZZ = t(Z)%*%Z
    bt.est = solve(ZZ, t(Z) %*%Y)
 }
 if (method=='LSP') {
   G = diag(ncol(Z))
   H = rep(0,ncol(Z))
   bt.est = lsei(A=Z,B=Y,G=G, H=H, verbose=FALSE)$X
 }
 if (method=='EM'){
    csum = colSums(Z)  ## should be =1, but here just in case
    for (i in 1:maxiter){
      X1 = t(t(Z)*beta0)                ## expected per cell under current est
       rsum = rowSums(X1)
      X2 = X1/ifelse(rsum>0, rsum,1)   ## normalize each row to get probabilities; avoid div by zero
     beta = c(t(X2) %*% Y) /csum        ## update:add up info across rows, then normalize by colSum
     err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
     if (max(err) < maxerr) break
     beta0 = beta
    }  ## end for i
  bt.est = beta
  } ## endif method=EM
 return(matrix(bt.est, nrow=n, byrow=TRUE))
}

## .. given (Tx expression, data Ymat, sample size n=ncol(Ymat)): estimate X
## .............. OLS: note Z is potentially a big matrix!!
##
estim.X = function(BETA,Ymat,method='LSP'){
 Y = c(Ymat)
 n = ncol(Ymat)    ## also = nrow(BETA)
 q = nrow(Ymat)    ## number of eqClasses per sample
 p = ncol(BETA)   ## number of transcripts
 sdy= sd(Y)
   Y1= Y/sdy       ## standardize to make lsp work better!!
 Z = kronecker(BETA, diag(q))
   Z = Z/sdy  ## normalize size
 if (method=='OLS'){
   ZZ = t(Z)%*%Z
   X = solve(ZZ, t(Z) %*%Y1)
 }
 if (method=='LSP'){
  G = diag(ncol(Z))
  H = rep(0,ncol(Z))
  X = lsei(A=Z,B=Y1,G=G, H=H, verbose=FALSE)$X
 }
 # normalize
 tmp = matrix(X, ncol=p)    
   X = t(t(tmp)/colSums(tmp))
 return(X)
}

## ................  EM based
## 
## calculate without using big matrix Z, but use for-loop
## the background pattern X.em is in matrix form
## .. given BETA= matrix of Tx expression, ** samples in rows ***
## Y = data vector, sample by sample
## n  = sample size n 
## ==> update X
#= X.em= matrix of starting values with OLS+
estim.X.em0 = function(X.em,BETA,Ymat){
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)  ## sample size
  p = ncol(BETA)  ## number of transcripts
 Ymat = matrix(Y,ncol=n)
 Ycell = NULL
 for (j in 1:n){
  bt = BETA[j,]                    ## expression-vector of j'th sample
  y = Ymat[,j]                   ## data vector from j'th sample
  x1 = t(t(X.em)*bt)             ## expected per cell
    x2 = x1/rowSums(x1)          ## normalize each row to get probabilities
    ycell = c(x2* y) 
    Ycell = rbind(Ycell, ycell)
  }
  BTsum = rep(colSums(BETA), rep(q,p))
  X.em = matrix(colSums(Ycell)/BTsum, ncol=p)
    X.em = t(t(X.em)/colSums(X.em))  ## normalize to sum to 1
  # print(round(c(X.em),2))
  return(X.em)
}

##
## EM without for-loop over sample size n, only over p
##
estim.X.em = function(X.em,BETA,Ymat){
  X0 = X.em  ## starting value
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X.em)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)     ## sample size
  p = ncol(BETA)  ## number of transcripts
 BTsum = rep(colSums(BETA), rep(q,p)) ## NOTE: BTsum might be zero!!
 x1 = NULL
 for (j in 1:p) x1=cbind(x1,kronecker(BETA[,j],X.em[,j]))
 sumx1 = rowSums(x1) 
   sumx1.pos= ifelse(sumx1>0, sumx1, 0.1)
 x2 = x1/sumx1.pos
 Ycell= x2*Y
 tmp = NULL 
 for (j in 1:p){
   tmp = c(tmp,colSums(matrix(Ycell[,j], nrow=n,byrow=TRUE)))
  }
 ## NOTE BTsum might be zero!!
 ## if csum=0, return to the original X0
 X.em = matrix(tmp/ifelse(BTsum>0, BTsum, 0.1), ncol=p)
   csum = colSums(X.em)
   X.em = t(t(X.em)/ifelse(csum>0,csum,0.1))  ## normalize to sum to 1
 for (j in 1:length(csum)){if (csum[j]< 0.00001) X.em[,j]=X0[,j]}
 return(X.em)
}

getSE <-function(X,b,y){
  #general for a matrix of y
    yhat=X %*% t(b)
    err=y-yhat
    W = 1/(yhat+1)  ## weight matrix
    SE = NULL
    for (j in 1:ncol(W)){
      WX = (W[,j]*X)
      XWX = t(X) %*% WX
      XWXi = solve(XWX)
      XWeWX = t(WX) %*% (err[,j]^2 *WX)
      mat = XWXi %*% XWeWX %*% XWXi
      se = sqrt(diag(mat))
      SE = rbind(SE,se)
    }
    return(SE)
}



estim_chunk <- function(Y,DM0,maxiter.X=1000,maxerr=0.01,lim=0.01,maxiter.bt=1000,modify=TRUE){
#  maxiter.X=1000
#  maxerr=0.01
#  lim=0.01
#  maxiter.bt=1000
#  modify=TRUE
  cos_sim=function(a,b) return(sum(a*b)/sqrt(sum(a^2)*sum(b^2)))

  X.y= Y

  DM0_ncol=sapply(DM0,ncol)
  DM0_nrow=sapply(DM0,nrow)

  poorCondX=NULL
  for (m in 1:length(DM0))
    if (DM0_ncol[m]==2)
  {
    X0=DM0[[m]]
    isPoor=TRUE
    for (i in 1:ncol(X0)){
      sim=apply(X0,2,cos_sim,X0[,i])
      if (sum(sim < 0.99)){ 
        isPoor=FALSE
        break()          
      }
    }
    if (isPoor) poorCondX=c(poorCondX,m)
  }

  adjust_esp=2;
  final <- NULL
  for(m in 1:length(Y)){
    Ymat  = X.y[[m]]
    X0=DM0[[m]]

    if(nrow(X0) == 1){
#      final[[m]] <- data.frame(Samples = colnames(Ymat), Count = Ymat[,colnames(Ymat)])
      rownames(Ymat)=colnames(X0)
      Ymat=t(Ymat)
      final=cbind(final,Ymat)
      
    } else {

      #CRP with only 2 txs
      needAdjustX=FALSE
      if (m %in% poorCondX){
        rs=rowSums(Ymat)
        rs=rs/sum(rs)
        sim=apply(X0,2,cos_sim,rs)
        e_sim=max(sim)-min(sim)
        e_sim
        bestTx=which.max(sim)
#        adjID=which.min(X0[,bestTx])
        adjID=which(X0[,bestTx]>0) #these eqclasses are supposed to be expressed
#        adjID=which(X0[,bestTx]>0 | rs>0)
        if (length(adjID)>0){
        #if (sum(rs[adjID]==0)>0){          
#          if (X0[adjID,bestTx] > max(X0[adjID,-bestTx])){
            #add a small value to Y
            needAdjustX=TRUE
            Ymat[adjID,]=Ymat[adjID,]+adjust_esp
  #          Ymat=Ymat+adjust_esp
#          }
        }
      }

      ### estimation
      X.est0 = X0  # starting design matrix
      Ysum= colSums(Ymat)
      beta0 = rep(Ysum/ncol(X0), rep(ncol(X0), length(Ysum))) # starting espression
      ## alternating estimation
      for (i in 1:maxiter.X){
        BT = estim.BETA(X.est0, Ymat, beta0=beta0, maxiter=maxiter.bt, modify=modify);
        X.est = estim.X.em(X.est0,BT,Ymat); ## print(rbind(BT,X.est))
        ## check convergence
        err = abs(X.est-X.est0)/ifelse(abs(X.est0)>lim, X.est0, lim)  ## X.est0 could be zero
        if (max(err) < maxerr) break
        ## update for next iterate
        X.est0 = X.est
        beta0=c(t(BT))   ## NOTE: BT is a matrix, saved byrow=T
      }
      colnames(BT)=colnames(X0)
      ### estimation - done
      if (needAdjustX){
        #need to subtract to the esp value
        #Ymat[adjID,]=Ymat[adjID,]+adjust_esp
        deductVal=length(adjID)*adjust_esp
        BT2=apply(BT,1,function(x){
          d=deductVal
          p=x/sum(x)
          x=x-d*p
          return(x)
        })
        BT=t(BT2)
#        yhat=X0 %*% t(BT)
#        err=y-yhat
      }

      final=cbind(final,BT)
#      final[[m]] <- data.frame(Samples = colnames(Ymat), Count = BT)
    }
  }

  rownames(final)=colnames(X.y[[1]])

  cat("\n estimate a chunk: Done!")
  return(final)
}



estim_chunk2_1 <- function(flist,DM0,maxiter.X=1000,maxerr=0.01,lim=0.01,maxiter.bt=1000,modify=TRUE){
#  maxiter.X=1000
#  maxerr=0.01
#  lim=0.01
#  maxiter.bt=1000
#  modify=TRUE
  Yfull=list()
  for (i in 1:length(flist)){
    cat(" ",i)
    load(flist[i])
    Yfull[[i]]=Y
  }

  cos_sim=function(a,b) return(sum(a*b)/sqrt(sum(a^2)*sum(b^2)))

  DM0_ncol=sapply(DM0,ncol)
  DM0_nrow=sapply(DM0,nrow)

  poorCondX=NULL
  for (m in 1:length(DM0))
    if (DM0_ncol[m]==2)
  {
    X0=DM0[[m]]
    isPoor=TRUE
    for (i in 1:ncol(X0)){
      sim=apply(X0,2,cos_sim,X0[,i])
      if (sum(sim < 0.99)){ 
        isPoor=FALSE
        break()          
      }
    }
    if (isPoor) poorCondX=c(poorCondX,m)
  }

  adjust_esp=2;
  final <- NULL
  for(m in 1:length(DM0)){
    cat(" ",m)
    X0=DM0[[m]]
    Ymat=NULL
    for (i in 1:length(flist)){      
      #load(flist[i])
      Y=Yfull[[i]]
      Ymat = cbind(Ymat,Y[[m]])
    }

    ls_crp=colSums(Ymat)
    if (sum(ls_crp)>0){

    if(nrow(X0) == 1){
      rownames(Ymat)=colnames(X0)
      Ymat=t(Ymat)
      final=cbind(final,Ymat)
      
    } else {

      allCellID=colnames(Ymat)
      pick=ls_crp > 0
      table(pick)
      Ymat=Ymat[,pick,drop=FALSE]

      #CRP with only 2 txs
      needAdjustX=FALSE
      if (m %in% poorCondX){
        rs=rowSums(Ymat)
        rs=rs/sum(rs)
        sim=apply(X0,2,cos_sim,rs)
        e_sim=max(sim)-min(sim)
        e_sim
        bestTx=which.max(sim)
#        adjID=which.min(X0[,bestTx])
        adjID=which(X0[,bestTx]>0) #these eqclasses are supposed to be expressed
#        adjID=which(X0[,bestTx]>0 | rs>0)
        if (length(adjID)>0){
        #if (sum(rs[adjID]==0)>0){          
#          if (X0[adjID,bestTx] > max(X0[adjID,-bestTx])){
            #add a small value to Y
            needAdjustX=TRUE
            Ymat[adjID,]=Ymat[adjID,]+adjust_esp
  #          Ymat=Ymat+adjust_esp
#          }
        }
      }

      ### estimation
      X.est0 = X0  # starting design matrix
      Ysum= colSums(Ymat)
      beta0 = rep(Ysum/ncol(X0), rep(ncol(X0), length(Ysum))) # starting espression
      ## alternating estimation
      for (i in 1:maxiter.X){
        BT = estim.BETA(X.est0, Ymat, beta0=beta0, maxiter=maxiter.bt, modify=modify);
        X.est = estim.X.em(X.est0,BT,Ymat); ## print(rbind(BT,X.est))
        ## check convergence
        err = abs(X.est-X.est0)/ifelse(abs(X.est0)>lim, X.est0, lim)  ## X.est0 could be zero
        if (max(err) < maxerr) break
        ## update for next iterate
        X.est0 = X.est
        beta0=c(t(BT))   ## NOTE: BT is a matrix, saved byrow=T
      }
      colnames(BT)=colnames(X0)
      ### estimation - done
      if (needAdjustX){
        #need to subtract to the esp value
        #Ymat[adjID,]=Ymat[adjID,]+adjust_esp
        deductVal=length(adjID)*adjust_esp
        BT2=apply(BT,1,function(x){
          d=deductVal
          p=x/sum(x)
          x=x-d*p
          return(x)
        })
        BT=t(BT2)
#        yhat=X0 %*% t(BT)
#        err=y-yhat
      }
      rownames(BT)=colnames(Ymat)
      BT3=matrix(0,nrow=length(allCellID),ncol=ncol(BT))
      rownames(BT3)=allCellID
      colnames(BT3)=colnames(BT)
      BT3[rownames(BT),]=BT
      final=cbind(final,BT3)
    }
   }
  }


  cat("\n estimate a chunk: Done!")
  return(final)
}


estim_chunk2 <- function(flist,DM0,maxiter.X=1000,maxerr=0.01,lim=0.01,maxiter.bt=1000,modify=TRUE){
#  maxiter.X=1000
#  maxerr=0.01
#  lim=0.01
#  maxiter.bt=1000
#  modify=TRUE
  packages <- c("data.table")
  
  if(!all(packages%in% rownames(installed.packages()))){
    needed <- packages[which(!packages %in% rownames(installed.packages()))]
    lapply(needed, install.packages, repos = "https://cran.r-project.org")
  }
  library("data.table")

  Yfull=list()
  for (i in 1:length(flist)){
    cat(" ",i)
    load(flist[i])
    Yfull[[i]]=Y
  }

  cos_sim=function(a,b) return(sum(a*b)/sqrt(sum(a^2)*sum(b^2)))

  DM0_ncol=sapply(DM0,ncol)
  DM0_nrow=sapply(DM0,nrow)

  poorCondX=NULL
  for (m in 1:length(DM0))
    if (DM0_ncol[m]==2)
  {
    X0=DM0[[m]]
    isPoor=TRUE
    for (i in 1:ncol(X0)){
      sim=apply(X0,2,cos_sim,X0[,i])
      if (sum(sim < 0.99)){ 
        isPoor=FALSE
        break()          
      }
    }
    if (isPoor) poorCondX=c(poorCondX,m)
  }

  adjust_esp=2;
#  final <- NULL
  allCellID=NULL
  for (i in 1:length(flist)){      
    allCellID = c(allCellID,colnames(Yfull[[i]][[1]]))
  }

#  allTx=NULL
#  for(m in 1:length(DM0)){
#    allTx=c(allTx,colnames(DM0[[m]]))
#  }

  final=data.table()
#  final=matrix(0,nrow=ncol(Yfull[[1]]),ncol=length(allTx))
#  rownames(final)=allCellID
#  colnames(final)=allTx
#  for(m in 1:length(DM0)){
  final=foreach(m = 1:length(DM0),.combine=rbind) %dopar% {
    cat(" ",m)
    X0=DM0[[m]]
    Ymat=NULL
    for (i in 1:length(flist)){      
      #load(flist[i])
      Y=Yfull[[i]]
      Ymat = cbind(Ymat,Y[[m]])
    }

    ls_crp=colSums(Ymat)
    if (sum(ls_crp)>0){

    if(nrow(X0) == 1){
      rownames(Ymat)=colnames(X0)
      #Ymat=t(Ymat)
      #final=cbind(final,Ymat)
      Ymat=data.table(Ymat)
      Ymat$Transcript=rownames(Ymat)
      #final=rbind(final,Ymat)
      return(Ymat)
    } else {      
      pick=ls_crp > 0
      table(pick)
      Ymat=Ymat[,pick,drop=FALSE]

      #CRP with only 2 txs
      needAdjustX=FALSE
      if (m %in% poorCondX){
        rs=rowSums(Ymat)
        rs=rs/sum(rs)
        sim=apply(X0,2,cos_sim,rs)
        e_sim=max(sim)-min(sim)
        e_sim
        bestTx=which.max(sim)
#        adjID=which.min(X0[,bestTx])
        adjID=which(X0[,bestTx]>0) #these eqclasses are supposed to be expressed
#        adjID=which(X0[,bestTx]>0 | rs>0)
        if (length(adjID)>0){
        #if (sum(rs[adjID]==0)>0){          
#          if (X0[adjID,bestTx] > max(X0[adjID,-bestTx])){
            #add a small value to Y
            needAdjustX=TRUE
            Ymat[adjID,]=Ymat[adjID,]+adjust_esp
  #          Ymat=Ymat+adjust_esp
#          }
        }
      }

      ### estimation
      X.est0 = X0  # starting design matrix
      Ysum= colSums(Ymat)
      beta0 = rep(Ysum/ncol(X0), rep(ncol(X0), length(Ysum))) # starting espression
      ## alternating estimation
      for (i in 1:maxiter.X){
        BT = estim.BETA(X.est0, Ymat, beta0=beta0, maxiter=maxiter.bt, modify=modify);
        X.est = estim.X.em(X.est0,BT,Ymat); ## print(rbind(BT,X.est))
        ## check convergence
        err = abs(X.est-X.est0)/ifelse(abs(X.est0)>lim, X.est0, lim)  ## X.est0 could be zero
        if (max(err) < maxerr) break
        ## update for next iterate
        X.est0 = X.est
        beta0=c(t(BT))   ## NOTE: BT is a matrix, saved byrow=T
      }
      colnames(BT)=colnames(X0)
      ### estimation - done
      if (needAdjustX){
        #need to subtract to the esp value
        #Ymat[adjID,]=Ymat[adjID,]+adjust_esp
        deductVal=length(adjID)*adjust_esp
        BT2=apply(BT,1,function(x){
          d=deductVal
          p=x/sum(x)
          x=x-d*p
          return(x)
        })
        BT=t(BT2)
#        yhat=X0 %*% t(BT)
#        err=y-yhat
      }
      rownames(BT)=colnames(Ymat)
      BT3=matrix(0,nrow=length(allCellID),ncol=ncol(BT))
      rownames(BT3)=allCellID
      colnames(BT3)=colnames(BT)
      BT3[rownames(BT),]=BT
      #final=cbind(final,BT3)
      BT3=t(BT3)
      BT3=data.table(BT3)
      BT3$Transcript=colnames(BT)
#      Ymat=data.table(Ymat)
#      Ymat$Transcript=rownames(Ymat)
#      final=rbind(final,BT3)
      return(BT3)
    }
   } #if (sum(ls_crp)>0){
    return (NULL)
  }


  cat("\n estimate a chunk: Done!")
  return(final)
}

