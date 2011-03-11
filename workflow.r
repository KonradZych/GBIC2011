#copy-paste into R for analysis
workflow.appriopriateMarkers <- function(){
	setwd("D:/data")
	library(basicQtl)
	library(qtl)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix,margin=0.5,genotypes=c(1,0),overlapInd=0, verb=T)
	brassicaFlipped <- flipMarkers(brassica_genotypes)
	cross <- orderedCross(brassica_genotypes,verbose=T)
	plot.rf(formLinkageGroups(cross,reorgMarkers=F))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix,margin=0.5, overlapInd=0, verb=T)
	brassica_genotypes <- switchMatrixValues(brassica_genotypes,before=c("A","B"),after=c(0,1))
	brassicaReco <- recombinationCount(brassica_genotypes)
	brassicaRecoflipped <- recombinationCount(brassica_genotypes,flip=1)
}