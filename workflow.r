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
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2)
	crossSimulated <- toGenotypes(ril,use="simulated",verbose=TRUE,debugMode=1,treshold=0.05,margin=15)
	
	setwd("D:/data/")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2)
	crossMean <- toGenotypes(ril,use="simulated",treshold=0.01,overlapInd = 0, proportion = 50, margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	c2 <- crossMean
	plot.map(cleanMap(c2),5)
	crossMean <- toGenotypes(ril,use="simulated",treshold=0.01,overlapInd = 0, proportion = 75, margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	switchChromosomes.internal
	crossMap <- toGenotypes(ril,use="map",treshold=0.01,overlapInd = 0, proportion = 50, margin=15, verbose=TRUE, debugMode=1)
	
	
	setwd("D:/data/bc")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2)
	crossMean <- toGenotypes(ril,use="simulated",treshold=0.01,overlapInd = 0, proportion = 75, margin=15, minChrLength=5, verbose=TRUE, debugMode=1)
	
	png("tolerateordereddefault.png", width=1200, height=1200)
	plot.rf(formLinkageGroups(cross,reorgMarkers=TRUE))
	dev.off()
 }

 sum(names(ril$parental$shrinkt[which(ril$parental$shrinkt < 0.01)]) %in% rownames(ril$parental$phenotypes[which(ril$parental$RP$pval[1]<0.01),]))
 
 sum(ril$parental$shrinky[which(ril$parental$shrinky<0.01)] %in% ril$parental$shrinkt[which(ril$parental$shrinkt<0.01)])
 
 upParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < 0.01),]
 downParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] < 0.01),]
 ndownParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] > 0.5),]
 nupParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] > 0.5),]
pare <- rbind(upParental[1:50,],downParental[1:50,],ndownParental[1:100,],nupParental[1:100,])
rils <- ril$rils$phenotypes[rownames(pare),]
 write.table(rils,file="children_phenotype.txt",sep="\t")
 write.table(pare,file="parental_phenotype.txt",sep="\t")
 m<-ril$rils$map[which(ril$rils$map[,3] %in% rownames(ril$rils$phenotypes[1:10,])),]
 write.table(m,file="children_mape.gff",sep="\t")
 
datasy<-read.csv(skip=200, sep="\t", fill=TRUE, file="GDS1115_full.soft")
 
for(i in rownames(m)){
	if(n%%100==0) print(n)
	m[i,1] <- sum(file2[,1]==i)
	n<-n+1
}

6157  markers
132 individuals
 
 
cleaner <- function(ril){
	for(i in rownames(ril$parental$phenotypes)){
		if()
	}
}



data <- read.table(file, sep=",",colClasses="character", fill=TRUE, blank.lines.skip=TRUE, comment.char="")

yoyo <- function(cross){
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