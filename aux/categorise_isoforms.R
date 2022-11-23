## 24Mar2021/Nghia: create annotation and mapping of paralogs, nonparalogs, singleton
# - input: sqlite and xmatrix
# - ouput: paralogs, nonparalogs, singleton mapping
# - command: Rscript categorise_isoforms.R sqlite=hg38.refGene.sqlite xmatrix=Xmatrix.RData out=isoform_groups_UCSC_hg38_mapper.RData

#24Mar2021/Nghia: fix the bug to make sure no genes are not counted twice in genes.tx.map.all.final
#11Mar2021/Nghia: fix the bug that some isoforms having no corespoding genes, so use the names of the isoforms instead
#18Mar2021:fix the bug that replaces gene names of non-paralog and singleton by paralog gene names


args = commandArgs(trailingOnly=TRUE)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="sqlite") sqliteFn=res[2]
  if (res[1]=="xmatrix") xmatrixFn=res[2]
  if (res[1]=="out") outFn=res[2]
}

library("GenomicFeatures")

anntxdb <- loadDb(sqliteFn)
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME"), keytype = "GENEID")))
genes.tx.map=genes.tx.all$GENEID
names(genes.tx.map)=genes.tx.all$TXNAME

genes.tx.all.export=genes.tx.all[,c("TXNAME","GENEID","GENEID")]
tx.all.ID=unique(genes.tx.all$TXNAME)

###################
################ categorise isoforms and genes
load(xmatrixFn)

c2t=sapply(names(CRP), strsplit, split=' ') 
c_len=lengths(c2t)
#map from tx to crp
t2c=unlist(c2t,use.names=FALSE)
names(t2c)=rep(names(c2t), lengths(c2t))
t2c=tapply(names(t2c),t2c,c)
#summary(lengths(t2c))
pick=which(lengths(t2c)>1)
#length(pick)

#singleton: genes or isoforms from genes with a single isoform
#gene with single isoform - singleton
singleton_gene=table(genes.tx.all$GENEID)
singleton_gene=singleton_gene[singleton_gene==1]
singleton_gene=names(singleton_gene)
singleton_isoform=genes.tx.all$TXNAME[genes.tx.all$GENEID%in%singleton_gene]

#non-paralogs: isoforms from genes with multiple isoforms and not merged
nonpara_isoform=unlist(lapply(CCRP, colnames))
names(nonpara_isoform)=NULL
pick=nonpara_isoform %in% genes.tx.all$TXNAME
nonpara_isoform=nonpara_isoform[pick]
pick=nonpara_isoform %in% singleton_isoform
nonpara_isoform=nonpara_isoform[!pick]
#length(nonpara_isoform)
nonpara_gene=genes.tx.all$GENEID[genes.tx.all$TXNAME%in%nonpara_isoform]

#isoform paralogs: isoforms from genes with multiple isoforms and merged with other isoforms
para_isoform=unlist(lapply(CCRP, colnames))
names(para_isoform)=NULL
pick=para_isoform %in% genes.tx.all$TXNAME
para_isoform=para_isoform[!pick]
#length(para_isoform)

# gene paralogs: genes from paralog isoforms with merged with other genes
tx_sep=strsplit(para_isoform," ")
para_gene.list=lapply(tx_sep,function(x){
  sapply(x,function(y) {
    z=genes.tx.map[y]
    if (is.na(z)) z=y
    return(z)
  })
})
para_gene.list=lapply(para_gene.list,unique)
para_gene.list=lapply(para_gene.list,sort)
para_gene.len=lengths(para_gene.list)

#genes from isoform paralogs
isopara_gene.names=sapply(para_gene.list,paste,collapse=" ")
names(isopara_gene.names)=para_isoform
#update the mapping
genes.tx.map.all=c(genes.tx.map,isopara_gene.names)
#head(isopara_gene.names)

#some genes belong to more than 1 gene paralog - merge them together
myclusts=1:length(para_gene.list)
for (i in 1:(length(para_gene.list)-1))
  if (para_gene.len[i]>1){
  g1=para_gene.list[[i]]
  for (j in (i+1):length(para_gene.list))
    if (para_gene.len[j]>1){
      g2=para_gene.list[[j]]
      if (length(intersect(g1,g2))>0) myclusts[j]=myclusts[i]
    }
}
x=table(myclusts)
x=as.integer(names(which(x>1)))
for (i in x){
  pick=which(myclusts==i)
  m=sort(unique(unlist(para_gene.list[pick])))
  for (j in pick) para_gene.list[[j]]=m
}
para_gene=sapply(para_gene.list,paste,collapse=" ")
names(para_gene)=para_isoform

### there are some singletons belong to paralogs, need to exclude them from the singleton sets
sharedGenes=intersect(singleton_gene,unlist(para_gene.list))
singleton_gene=setdiff(singleton_gene,sharedGenes)
tx_sep=strsplit(para_isoform," ")
sharedTxs=intersect(singleton_isoform,unlist(tx_sep))
singleton_isoform=setdiff(singleton_isoform,sharedTxs)


#consider the isolate genes
para_gene.names=unique(unlist(para_gene.list[para_gene.len>1]))
genes.tx.map2=genes.tx.map
#18Mar2021:fix the bug that replaces gene names of non-paralog and singleton by paralog gene names
#pick=which(genes.tx.map2 %in% para_gene.names)
pick=which(genes.tx.map2 %in% para_gene.names & !(names(genes.tx.map2) %in% c(nonpara_isoform,singleton_isoform)))
for (i in pick){
  g=genes.tx.map2[i]
  j=grep(g,para_gene)
  genes.tx.map2[i]=para_gene[j[1]]
}
genes.tx.map.all.final=c(genes.tx.map2,para_gene)

#### process paralog genes - when a paralog from 1 gene can belong to a paralog from multiple genes. Remeber the cases of multiple genes have been processed before
genes.tx.map.raw=genes.tx.map.all.final

genes=genes.tx.map.all.final[c(nonpara_isoform,para_isoform)]
gene_set=strsplit(genes," ")
gene_set.len=lengths(gene_set)
g2 = which(gene_set.len>1)
g1 = which(gene_set.len==1)
myid=NULL
for (i in 1:length(g2)){
  p1=g2[i]
  mygeneset=gene_set[[p1]]
  pick=which(genes[g1] %in% mygeneset)
  if (length(pick)>0){
    paraList=names(genes[g1][pick])
    genes.tx.map.all.final[paraList]=genes[p1]#update
    myid=c(myid,i)
  }
}



save(para_isoform,nonpara_isoform,singleton_isoform,genes.tx.map.all.final,genes.tx.map.raw, file=outFn)
