convertToGenotypes <- function(ril, overlapInd, treshold, proportion, margin, verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	if(length(proportion)==2){
		genotypesUp <- c(0,1)
		genotypesDown <- c(1,0)
	}else if(length(proportion)==3){
		genotypesUp <- c(0,1,2)
		genotypesDown <- c(2,1,0)
	}
	m <- NULL
	upParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < treshold),]
	downParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] < treshold),]
	upRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(upParental)),]
	downRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(downParental)),]
	for(x in rownames(upRils)){
		m <- rbind(m,splitRow.internal(x, upRils, upParental, overlapInd, proportion, margin, ril$parental$groups, genotypesUp))
	}
	for(x in rownames(downRils)){
		m <- rbind(m,splitRow.internal(x, downRils, downParental, overlapInd, proportion, margin, ril$parental$groups, genotypesDown))
	}
	ril$rils$genotypes$simulated <- m
	colnames(ril$rils$genotypes$simulated) <- colnames(upRils)
	rownames(ril$rils$genotypes$simulated) <- c(rownames(upRils),rownames(downRils))
	invisible(ril)
}

splitRow.internal <- function(x, rils, parental, overlapInd, proportion, margin, groups, genotypes){
	result <- rep(0,length(rils[x,]))
	A <- parental[which(rownames(parental) == x),which(groupLabels==0)]
	B <- parental[which(rownames(parental) == x),which(groupLabels==1)]
	splitVal <- mean(mean(A,na.rm=TRUE),mean(B,na.rm=TRUE))
	if(length(proportion)==2){
		result[which(rils[x,] > splitVal)] <- genotypes[1]
		result[which(rils[x,] < splitVal)] <- genotypes[2]
		result[which(rils[x,] == splitVal)] <- NA
		### this should be merged!
		result <- filterRow.internal(result)
	}else if(length(proportion)==3){
		
	}
	invisible(result)
}
