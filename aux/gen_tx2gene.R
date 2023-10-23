## Generate the data of the mapping tx to gene for GENCODE, ENSEMBL and other annotations
# - input: cdna file, gtf file, clean (optional), and anntype (optional)
# - ouput: sqlite, txp2gene and cdnaout which contain only common transcripts between gtf and cdna
# - command: Rscript gen_tx2gene.R cdna=Homo_sapiens.GRCh38.cdna.all.fa gtf=Homo_sapiens.GRCh38.106.gtf anntype=ENSEMBL sqlite=Homo_sapiens.GRCh38.106.sqlite cdnaout=Homo_sapiens.GRCh38.106.cdna.remain.fa out=txp2gene_GRCh38.106.tsv

### Nghia/12Sep2023: 
# - revise to process GENCODE annotation

#### Nghia/23Oct2023:
# - add an extra parameter: clean which allows removing the transcripts from cdna but not found in the gtf (clean=1, default) or not (clean=0). We observe that those transcripts without corresponding genes are mainly non-coding transcripts (NR_* in Refseq), if they are your interesting transcripts, please use clean=0.



anntype="ENSEMBL"
clean=1
args = commandArgs(trailingOnly=TRUE)
for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="cdna") cdnaFn=res[2]
	if (res[1]=="anntype") anntype=res[2]
	if (res[1]=="clean") clean=as.integer(res[2])
	if (res[1]=="gtf") gtfFn=res[2]
	if (res[1]=="sqlite") sqliteFn=res[2]
	if (res[1]=="cdnaout") cdnaoutFn=res[2]
	if (res[1]=="out") outFn=res[2]
}

cat("\n------------------------------")
cat("\n Running with the following parameter settings: ")
	cat("\n cdna=",cdnaFn)
	cat("\n anntype (ENSEMBL (default), GENCODE, or OTHER)=",anntype)
	cat("\n clean (1: default, remove transcripts not in gtf; or 0: include all transcripts )=",clean)
	cat("\n gtf=",gtfFn)
	cat("\n sqlite=",sqliteFn)
	cat("\n cdnaout=",cdnaoutFn)
	cat("\n out=",outFn)
cat("\n------------------------------")

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
if (anntype=="ENSEMBL") cdna_tx=sapply(names(fasta_tx),function(x) unlist(strsplit(x," "))[1])
if (anntype=="GENCODE") cdna_tx=sapply(names(fasta_tx),function(x) unlist(strsplit(x,"\\|"))[1])
names(cdna_tx)=NULL

if (anntype=="ENSEMBL")	cdna_tx=sapply(cdna_tx,function(x) unlist(strsplit(x,"\\."))[1]) #remove the version of the tx


if (clean==0){ #include all transcripts
	#get the concordant tx
	pick=genes.tx.all$TXNAME %in% cdna_tx  
	tx2gene=genes.tx.all[pick,]
	tx2gene=tx2gene[,c("TXNAME","GENEID")]

	p=cdna_tx %in% genes.tx.all$TXNAME
	cat("not in gtf: ",table(!p)) # some are not matched
	p1=which(!p)
	if (length(p1)>0){
		x=data.frame(TXNAME=cdna_tx[p1],GENEID=cdna_tx[p1])
		tx2gene=rbind(tx2gene,x)
	}


	tx.export.fasta=fasta_tx
	names(tx.export.fasta)=cdna_tx
	names(tx.export.fasta)=paste0(names(tx.export.fasta)," ",tx2gene$GENEID[match(names(tx.export.fasta),tx2gene$TXNAME)])
	writeXStringSet(tx.export.fasta, cdnaoutFn)

	
	write.table(tx2gene,file=outFn, col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")


}else{ #clean=1, remove transcripts if they are not in gtf
	#get the concordant tx
	pick=genes.tx.all$TXNAME %in% cdna_tx  
	tx2gene=genes.tx.all[pick,]
	tx2gene=tx2gene[,c("TXNAME","GENEID")]

	p=cdna_tx %in% genes.tx.all$TXNAME
	cat("not in gtf: ",table(!p)) # some are not matched

	tx.export.fasta=fasta_tx
	names(tx.export.fasta)=cdna_tx
	tx.export.fasta=tx.export.fasta[p]
	names(tx.export.fasta)=paste0(names(tx.export.fasta)," ",genes.tx.all$GENEID[match(names(tx.export.fasta),genes.tx.all$TXNAME)])
	writeXStringSet(tx.export.fasta, cdnaoutFn)

	tx.export.fasta=fasta_tx
	names(tx.export.fasta)=cdna_tx
	tx.export.fasta=tx.export.fasta[!p]
	names(tx.export.fasta)=paste0(names(tx.export.fasta)," ","Unknown-Gene-Name")
	writeXStringSet(tx.export.fasta, "tx_removed.fa")

	write.table(tx2gene,file=outFn, col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")
}



cat("Done!")
