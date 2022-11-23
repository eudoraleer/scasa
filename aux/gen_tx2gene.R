## nghiavtr/06July2022: generate mapping tx to gene
# - input: cdna and gtf file
# - ouput: sqlite, txp2gene and cdnaout which contain only common transcripts between gtf and cdna
# - command: Rscript gen_tx2gene.R cdna=Homo_sapiens.GRCh38.cdna.all.fa gtf=Homo_sapiens.GRCh38.106.gtf sqlite=Homo_sapiens.GRCh38.106.sqlite cdnaout=Homo_sapiens.GRCh38.106.cdna.remain.fa out=txp2gene_GRCh38.106.tsv
args = commandArgs(trailingOnly=TRUE)
for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="cdna") cdnaFn=res[2]
	if (res[1]=="gtf") gtfFn=res[2]
	if (res[1]=="sqlite") sqliteFn=res[2]
	if (res[1]=="cdnaout") cdnaoutFn=res[2]
	if (res[1]=="out") outFn=res[2]
}

suppressMessages(library("GenomicFeatures"))
suppressMessages(library("Biostrings"))

anntxdb <- makeTxDbFromGFF(file=gtfFn,
                 format="gtf",
                 dataSource=paste("Link to the source",sep=""),
                 organism="Homo sapiens")
saveDb(anntxdb,file=sqliteFn)

#anntxdb <- loadDb(sqliteFn)
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("TXNAME","GENEID"), keytype = "GENEID")))

#read fasta file
fasta_tx = readDNAStringSet(cdnaFn)
#get tx in cdna
cdna_tx=sapply(names(fasta_tx),function(x) unlist(strsplit(x," "))[1])
names(cdna_tx)=NULL
cdna_tx=sapply(cdna_tx,function(x) unlist(strsplit(x,"\\."))[1])
#get the concordant tx
pick=genes.tx.all$TXNAME %in% cdna_tx  
tx2gene=genes.tx.all[pick,]
tx2gene=tx2gene[,c("TXNAME","GENEID")]

p=cdna_tx %in% genes.tx.all$TXNAME
cat("not in gtf: ",table(p)) # some are not matched

tx.export.fasta=fasta_tx
names(tx.export.fasta)=cdna_tx
tx.export.fasta=tx.export.fasta[p]
names(tx.export.fasta)=paste0(names(tx.export.fasta)," ",genes.tx.all$GENEID[match(names(tx.export.fasta),genes.tx.all$TXNAME)])

writeXStringSet(tx.export.fasta, cdnaoutFn)


write.table(tx2gene,file=outFn, col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

cat("Done!")
