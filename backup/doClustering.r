fakeGenotypes <- function(len, outputFile, verbose=FALSE, debugMode=0){
	cat("faked,1,1,",file=outputFile,append=TRUE)
	cat(paste(sample(c(0,1),100,T),",",sep="",collapse=""),file=outputFile,append=TRUE)
	cat("\n",file=outputFile,append=TRUE)
}

fakePhenotypes <- function(len, outputFile, verbose=FALSE, debugMode=0){
	cat("faked,,,",file=outputFile,append=TRUE)
	cat(paste(sample(c(0.03,0.04,0.05),100,T),",",sep="",collapse=""),file=outputFile,append=TRUE)
	cat("\n",file=outputFile,append=TRUE)
}

doClustering <- function(genotypeMatrix, groups, iterations, outputFile, verbose=FALSE, debugMode=0){
	r <- bestClustering(genotypeMatrix,groups,iterations,verbose=verbose,debugMode=debugMode)
	r <- kmeans(r,groups)
	sorted <- sort(r[[1]],index.return=TRUE)
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose && debugMode==1) cat("writeGenotypes starting  for chromome",i,"out of",groups,".\n")
		writeGenotypes(genotypeMatrix[sorted[[2]][which(sorted[[1]]==i)],], i,outputFile, verbose, debugMode)
		el <- proc.time()
	}
	if(verbose && debugMode==2)cat("doClustering done in:",(el-sl)[3],"seconds.\n")
}