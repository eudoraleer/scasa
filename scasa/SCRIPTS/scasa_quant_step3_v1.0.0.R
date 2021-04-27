#!/usr/bin/Rscript
##############################################################################
#                  SCASA: QUANTIFICATION
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 3
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="workdir") workdir=as.character(res[2])
  if (res[1]=="txgroup") txgroup=as.character(res[2])
}

isoFn = paste(workdir,"/scasa_isoform_expression.RData",sep = "")
load(isoFn)
load(txgroup)

parasum = function(tx, count_isoform) {
     ncell = ncol(count_isoform)
     pick = tx %in% rownames(count_isoform)
     if (sum(pick)>1) {return(colSums(count_isoform[tx[pick],]))}
     if (sum(pick)==1){return(count_isoform[tx[pick],])}
     if (sum(pick)==0){return(rep(0, ncell))}
 }

scasa_genemat = function(scasa_iso, txmap)
{
  tx = rownames(scasa_iso)
  genes = txmap[tx]
  txlist = tapply(tx, genes, c)  # list of tx/paralogs
  ntx = sapply(txlist, length)
  tx1 = unlist(txlist[ntx==1])
  mat1 = scasa_iso[tx1,]
    rownames(mat1) = names(tx1)
  mat2 = sapply(txlist[ntx>1], parasum, count_isoform=scasa_iso)
    colnames(mat2) = names(txlist[ntx>1])
  gcount = rbind(mat1, t(mat2))
    colnames(gcount) = colnames(scasa_iso)
  return(gcount)
}

gene_count=scasa_genemat(isoform_count,genes.tx.map.all.final)

save(gene_count,file=paste(workdir,"/scasa_gene_expression.RData", sep = ""))
write.table(gene_count, paste(workdir,"/scasa_gene_expression.txt",sep = ""), quote = F, row.names = T, sep = "\t")
print("Estimate gene expression : done!")

