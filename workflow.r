#copy-paste into R for analysis
workflow.appriopriateMarkers <- function(){
	setwd("D:/data")
	library(basicQtl)
	library(qtl)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- toGenotypes(expressionMatrix,margin=0.5,genotypes=c(1,0),overlapInd=0, verb=T)
	recoMatrix <- recombinationCount(brassica_genotypes,flip=0)
	recoMatrixFlipped <- recombinationCount(brassica_genotypes,flip=1)
	res <- rep(0,length(rownames(recoMatrix)))
	names(res) <-  rownames(recoMatrix)
	for(x in 1:nrow(recoMatrix)){
		i <- recoMatrix[x,]
		j <- recoMatrixFlipped[x,]
		res[names(which(j[which(i>24)]<24))] <- res[names(which(j[which(i>24)]<24))]+1
	}
	hist(res)
	genotypes[which(res>780&&res<880),] <- 1-genotypes[which(res>780&&res<880),]
}

workflow.parental <- function(){
	setwd("D:/data/parental")
	library(pheno2geno)
	library(qtl)
	library(iqtl)
	expressionChildren <- as.matrix(read.table("Gene_quant.txt",sep=""))
	expressionParental <- as.matrix(read.table("Gene_parental.txt",sep=""))
	genotypesChildren <- as.matrix(read.table("Genotypes.txt",sep=" "))
	cross <- genotypesToCross(genotypesChildren,expressionChildren,verbose=TRUE,debugMode=2)
	correctionValues <-yoyo(cross)
	crosscor <- genotypesToCross(genotypesChildren,expressionCorrected ,verbose=TRUE,debugMode=2)
	correctionValues <-yoyo(crosscor)
	expressionCorrected <- expressionChildren - t(correctionValues[-147,])
	expressionChildrenCor <- cor(expressionChildren,use="pairwise.complete.obs")
	correctedChildrenCor <- cor(expressionCorrected,use="pairwise.complete.obs")
	
	parentalTtest <- -log10(t.test(expressionParental[,1:2],expressionParental[,3:4])[[3]])
}

readChildren <- function(){
	res1 <- as.matrix(read.table("Genotypes.txt",sep="\t",header=T))
	res1 <- t(res1)
	res2 <- res1[-1,]
	rownames(res2)<-rownames(res1[-1,], do.NULL = FALSE)
	colnames(res2)<-res1[1,]
	res2
}

yoyo <- function(cross){
  pheno <- pull.pheno(cross)
  cat("NROW " ,nrow(  pheno),"\n")
  cat("NCOL " ,ncol(  pheno),"\n")
  s <- proc.time()
  cormatrix <- NULL
  for(y in 1:ncol(pheno)){
    if(y %% 1000 == 0){
      e <- proc.time()
      cat("Done:",y,"/",ncol(pheno),"in",as.numeric(e[3]-s[3]),"secs\n")
      s <- e
    }
    oamean <- mean(pheno[,y])
    groupmeans <- lapply(lapply(batchlist, fun <- function(x){pheno[x,y]}),mean)
    diffmeans <- unlist(groupmeans) - oamean
    traitcorrection <- rep(0,nrow(pheno))
    for(x in 1:length(diffmeans)){
      traitcorrection[batchlist[[x]]] <- diffmeans[x]
    }
    cat(traitcorrection,"\n",sep="\t",file="tmpbatch.out",append=TRUE)
  }
  cormatrix <- t(read.table("tmpbatch.out",sep="\t"))
}