#!/usr/bin/Rscript
##############################################################################
#                  SCASA: GENERATE EQCLASS
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.0
#             Step: 1_1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2021-04-07
##############################################################################
packages <- c("plyr","data.table")

if(!all(packages%in% rownames(installed.packages()))){
  needed <- packages[which(!packages %in% rownames(installed.packages()))]
  lapply(needed, install.packages)
}
library("plyr")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
input_dir <- as.character(args[1])
output_dir <- as.character(args[2])
libsize_thres <- as.numeric(args[3])

if(!dir.exists(output_dir)) dir.create(output_dir)

cat("\n Export eqClasses of individual cells... ")
bus_output <- fread(paste(input_dir,"/output.correct.sort.bus.final.txt", sep = ""))
colnames(bus_output) <- c("Barcode","UMI_Code","eqClass","Count","Barcode_UMI")

libsize=bus_output[, .N, by = Barcode]
keepCells=libsize[which(libsize$N > libsize_thres),]
bus_output=bus_output[which(Barcode %in% keepCells$Barcode),]

eq <- fread(paste(input_dir,"/matrix.ec", sep = ""))
colnames(eq) <- c("eqClass","Transcript_Number")
eq$eqClass <- as.character(eq$eqClass)
eq$Transcript_Number <- as.character(eq$Transcript_Number)

ref <- read.table(paste(input_dir,"/transcripts.txt", sep = ""))
colnames(ref) <-c("Transcript_ID")
ref$Transcript_Number <- seq(0,(nrow(ref)-1),1)
#25Jan2021/Nghia: fix the bug when a total number of transcripts is greater than 100000(1e+05)
ref$Transcript_Number = as.character(as.integer(ref$Transcript_Number))

s <- strsplit(eq$Transcript_Number, split = ",")
eq <- data.frame(eqClass = rep(eq$eqClass, sapply(s, length)), Transcript_Number = unlist(s))
eq$Transcript_ID <- ref[match(eq$Transcript_Number,ref$Transcript_Number),"Transcript_ID"]

cell_umi <- bus_output[, .N, by = Barcode_UMI]
colnames(cell_umi)=c("Var1","Freq")

eq_size <- data.frame(table(eq$eqClass))
matchID=match(bus_output$eqClass,eq_size$Var1)
bus_output$eqSize <- eq_size[matchID,"Freq"]

#duplicated <- cell_umi[which(cell_umi$Freq > 1),]
matchID=match(bus_output$Barcode_UMI,cell_umi$Var1)
bus_output$Duplicate_Freq <- cell_umi[matchID,"Freq"]
bus_duplicated <- bus_output[which(bus_output$Duplicate_Freq > 1),]
bus_output <- bus_output[which(bus_output$Duplicate_Freq == 1),]
duplicated_unique <- unique(bus_duplicated$Barcode_UMI)

### 01Feb2021/Nghia: fast version to deal with duplicated Barcode_UMI
maxCountNum=bus_duplicated[,sum(Count==max(Count)),by=.(Barcode_UMI)]
## case1: only 1 eqclass with max count
case1_BU=maxCountNum[which(V1==1),]
case1=bus_duplicated[which(Barcode_UMI %in% case1_BU$Barcode_UMI),]
case1$id=seq(nrow(case1))
if (nrow(case1) > 0){
	pick=case1[,.(id[which(Count==max(Count))]),by=.(Barcode_UMI)]
	case1=case1[pick$V1,]
}

## case2:if there are more than 1 eclass with max count, then select the groups (Barcode_UMI) with only one row of max eqSize, otherwise (case3) excluded
case2_BU=maxCountNum[which(V1>1),]
case2=bus_duplicated[which(Barcode_UMI %in% case2_BU$Barcode_UMI),]
case2$id=seq(nrow(case2))
#keep ones with the max count (more than 1 in each group)
pick=case2[,.(id[which(Count==max(Count))]),by=.(Barcode_UMI)]
case2=case2[pick$V1,]
#then keep the groups (Barcode_UMI) with only one row of max eqSize
maxEqfreqNum=case2[,sum(eqSize==max(eqSize)),by=.(Barcode_UMI)]
pick=maxEqfreqNum[which(V1==1),]
case2=case2[which(Barcode_UMI %in% pick$Barcode_UMI),]
#finally get the row with the max eqSize
case2$id=seq(nrow(case2)) #update new id
if (nrow(case2)>0){
	pick=case2[,.(id[which(eqSize==max(eqSize))]),by=.(Barcode_UMI)]
	case2=case2[pick$V1,]
}

##add to bus_output
case1=case1[,c("Barcode","UMI_Code","eqClass","Count")]
case2=case2[,c("Barcode","UMI_Code","eqClass","Count")]
caseExtra=rbind(case1,case2)
case2=case1=NULL
bus_output=bus_output[,c("Barcode","UMI_Code","eqClass","Count")]
bus_output=rbind(bus_output,caseExtra)
caseExtra=NULL

bus_output$Count <- 1
selected_cells <- unique(bus_output$Barcode)

fun <-function(eqClass,Barcode){
	mycellID=Barcode[1]
	current_summary=data.frame(table(eqClass))
	colnames(current_summary)=c("eqClass","freq")
	pick=which(eq$eqClass %in% eqClass)
	current_final <- data.frame(Transcript = eq[pick,"Transcript_ID"],
	  Count = NA,
	  eqClass = eq[pick,"eqClass"])

	current_final$Count <- current_summary[match(current_final$eqClass,current_summary$eqClass),"freq"]
	current_final <- current_final[!is.na(current_final$Count),]
	current_final$Weight <- 1
	return(current_final)
}

cat("\n The number of cells: ",length(selected_cells))
sample_eqclass=bus_output[,fun(eqClass,Barcode),Barcode]
save(sample_eqclass,file=paste(output_dir,"/Sample_eqClass.RData", sep = ""))