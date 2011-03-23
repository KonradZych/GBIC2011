parentalSplit <- function(x,expressionChildren,expressionParental,verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1)if(x%%100==0)cat("parentalSplit starting withour errors in checkpoint for row:",x,"\n")
	s <- proc.time()
	expressionParentalRow <- expressionParental[which(rownames(expressionParental) %in% rownames(expressionChildren)[x]),]
	genotypeMatrixRow <- lapply(expressionChildren[x,],parentalSplitSub,expressionParentalRow)
	e <- proc.time()
	if(verbose && debugMode==2)if(x%%100==0)cat("parentalSplit for row:",x,"done in:",(e-s)[3],"seconds.\n")
	invisible(genotypeMatrixRow)
}

parentalSplitSub <- function(expressionChildrenElement,expressionParentalRow){
	distance1 <- abs(expressionChildrenElement-expressionParentalRow[1])
	distance2 <- abs(expressionChildrenElement-expressionParentalRow[2])
	if(distance1<=distance2){
		genotypeMatrixElement <- 0
	}else{
		genotypeMatrixElement <- 1
	}
	invisible(genotypeMatrixElement)
}

checkGeno <- function(genoMatrix1,genoMatrix2){
	genoMatrix1 <- mapMarkers(genoMatrix1,genoMatrix2,mapMode=1)
	genoMatrix1 <- mapMarkers(genoMatrix1,genoMatrix2,mapMode=2)
	genoMatrix2 <- mapMarkers(genoMatrix2,genoMatrix1,mapMode=1)
	genoMatrix2 <- mapMarkers(genoMatrix2,genoMatrix1,mapMode=2)
	print(dim(genoMatrix2))
	print(dim(genoMatrix1))
	correct <- sum(as.numeric(genoMatrix1)==as.numeric(genoMatrix2))/length(genoMatrix1)
	cat("Matrices are the same in:",correct,"%\n")
}