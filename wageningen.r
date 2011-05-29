############################################################################################################
#correctRowLoc.internal: probesLocation sub function, calculating mean possition of the gene and removes 
# begining and end information
# 
# genesRow - currently processed row
#
############################################################################################################
readBigFile.internal <- function(filename, rowBeg=NULL, rowEnd=NULL, )
{
	if(!file.exists(filename)) stop("File: ",filename,"doesn't exist.\n")
	if(is.null(rowBeg)){
		if(is.null(rowEnd)){
			result <- readAnnFile.internal(filename)
		}else{
			result <- readAnnFile.internal(filename,nrow=rowEnd)
		}
	}else{
		if(is.null(rowEnd)){
			result <- readAnnFile.internal(filename,skip=rowBeg-1)
		}else{
			result <- readAnnFile.internal(filename,skip=rowBeg-1,nrow=(rowEnd-rowBeg))
		}
	}
	invisible(result)
}

############################################################################################################
#correctRowLoc.internal: probesLocation sub function, calculating mean possition of the gene and removes 
# begining and end information
# 
# genesRow - currently processed row
#
############################################################################################################
readAnnFile.internal <- function(filename, ...){
	result <- read.table(filename,...)
	write.table(result[,c(3,4)],file="location_at.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,9,16,17,10,11,12)],file="parental_phenotypes.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,18,19,20,13,14,15)],file="parental_phenotypes_1.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,c(21:121)[which((21:121)%%2==1)])],file="children_phenotypes.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,c(21:121)[which((21:121)%%2==1)])],file="children_phenotypes_1.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
readAnnFile.internal <- function(filename, ...){
	result <- read.table(filename,...)
	write.table(result[,c(18,19,190,191)],file="parental_phenotypes_6h.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(20,21,192,193)],file="parental_phenotypes_ar.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(22,23,194,195)],file="parental_phenotypes_fresh.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(24,25,196,197)],file="parental_phenotypes_rp.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,5,26,35,36,37,52,56,58,63,69,70,75,84,85,88,92,94,96,108,109,112,113,114,115,128,134,142,147,148,149,155,157,160,161,165,166,167,168,173,175,178,180)],file="children_phenotypes_6h.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,5,31,33,34,39,43,47,50,53,55,60,65,68,76,77,78,79,82,91,93,99,101,107,116,124,126,129,130,137,139,144,146,152,153,164,169,172,176,179,184,187,188)],file="children_phenotypes_ar.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,5,27,28,32,38,40,41,45,46,51,59,61,64,66,72,74,80,86,87,89,95,98,102,105,117,118,120,121,123,132,133,136,141,143,150,159,170,171,174,181,182,189)],file="children_phenotypes_fresh.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	write.table(result[,c(4,5,29,30,42,44,48,49,54,57,62,67,71,73,81,83,90,97,100,103,104,106,110,111,119,122,125,127,131,135,138,140,145,151,154,156,158,162,163,177,183,185,186)],file="children_phenotypes_rp.txt",append=TRUE,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

readTomatoe <- function(){
	#cat("",file="parental_phenotypes_6h.txt",append=FALSE)
	#cat("",file="parental_phenotypes_ar.txt",append=FALSE)
	#cat("",file="parental_phenotypes_fresh.txt",append=FALSE)
	#cat("",file="parental_phenotypes_rp.txt",append=FALSE)
	#cat("",file="children_phenotypes_6h.txt",append=FALSE)
	#cat("",file="children_phenotypes_ar.txt",append=FALSE)
	#cat("",file="children_phenotypes_fresh.txt",append=FALSE)
	#cat("",file="children_phenotypes_rp.txt",append=FALSE)
	for(i in 20:100){
		skip <- (i-1)*10000+1
		nrows <- i*10000
		readAnnFile.internal("exp_ann_norm.txt",skip=skip,nrow=10000,sep="")
		cat(skip,nrows,"\n")
	}
}