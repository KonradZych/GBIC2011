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
	
	design <- sort(ril$rils$phenotypes[1,],decreasing=TRUE)
	weight <- c(rep(0,73),rep(0,73))
	modellikelihood(design[c(1:10,35000:35010),],rep(3,146),c(rep(0,73),rep(0,73)))
	
	setwd("D:/GenesForSeedQualityWageningen")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2,sep="\t")
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2,groupLabels=c(0,0,0,1,1,1))
	ril2 <- ril
	ril2$rils$phenotypes <- ril2$rils$phenotypes[,-1]
	ril2$rils$phenotypes <- matrix(as.numeric(ril2$rils$phenotypes),nrow(ril$rils$phenotypes),(ncol(ril$rils$phenotypes)-1))
	ril2$parental$phenotypes <- ril2$parental$phenotypes[,-1]
	ril2$parental$phenotypes <- matrix(as.numeric(ril2$parental$phenotypes),nrow(ril$parental$phenotypes),(ncol(ril$parental$phenotypes)-1))
	rownames(ril2$rils$phenotypes)<- 1:nrow(ril2$rils$phenotypes)
	colnames(ril2$rils$phenotypes)<- 1:ncol(ril2$rils$phenotypes)
	rownames(ril2$parental$phenotypes)<- 1:nrow(ril2$parental$phenotypes)
	colnames(ril2$parental$phenotypes)<- 1:ncol(ril2$parental$phenotypes)
	cross <- toGenotypes(ril2,use="simulated",treshold=0.05,overlapInd = 0, proportion = c(50,50), margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	
	setwd("D:/data/from slave/data")
	ril <- readFiles(verbose=TRUE,debugMode=2,sep="")
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2,groupLabels=c(0,0,1,1))
	cross <- toGenotypes(ril,use="simulated",treshold=0.05,overlapInd = 0, proportion = c(50,50), margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	
exp_p <- read.table(file="children_phenotypes_rp.txt" ,sep="\t")
exp_p <- read.table(file="parental_phenotypes_rp.txt" ,sep="\t")
	
	
	setwd("D:/data/bc")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,groupLabels = c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1),verbose=TRUE,debugMode=2)
	crossMean <- toGenotypes(ril,use="simulated",treshold=0.01,overlapInd = 0, proportion = c(25,50,25), minChrLength=5, verbose=TRUE, debugMode=1)
	
	png("tolerateordereddefault.png", width=1200, height=1200)
	plot.rf(formLinkageGroups(cross,reorgMarkers=TRUE))
	dev.off()
 }

ril <- convertToGenotypes.internal(ril, 0, 0.01, c(50,50), 15, TRUE, 1)







parental <- read.table(file="mapping_probes.txt",sep="\t",header=T,row.names=1)
children <- read.table(file="children.txt",sep="\t",header=T,row.names=1)
ril <- NULL
ril$rils$phenotypes <- children[,20:120]
ril$parental$phenotypes <- parental[,5:14]
ril <- preprocessData(ril, c(1,1,1,1,1,0,0,0,0,0))
cross <- toGenotypes(ril, use="simulated", splitMethod="EM", treshold=0.001, verbose=TRUE,debug=1)
