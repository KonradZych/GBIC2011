up <- c(0.3,0.55,0.3)
down <- c(0.2,0.45,0.2)
res <- NULL
for(i in 1:ncol(ril2$rils$EM)){
	cur <- ril2$rils$EM[c(1,2,3),i]
	cur2 <- ril2$rils$EM[c(1,3,2),i]
	cur3 <- ril2$rils$EM[c(2,1,3),i]
	printable <- FALSE
	if(all(cur<up)&&all(cur>down)){
		printable <- TRUE
	}else if(all(cur2<up)&&all(cur2>down)){
		printable <- TRUE
	}else if(all(cur3<up)&&all(cur3>down)){
		printable <- TRUE
	}
	if(printable){
		res <- c(res,i)
	}
	if(i%%100==0){
		print(ril2$rils$EM[,i])
		print(printable)
	}
}


convertToGenotypes.internal <- function(ril, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	if(length(proportion)==2){
		genotypesUp <- c(0,1)
		genotypesDown <- c(1,0)
	}else if(length(proportion)==3){
		genotypesUp <- c(0,1,2)
		genotypesDown <- c(2,1,0)
	}
	output <- NULL
	markerNames <- NULL 
	cat("--------1---------\n")
	upParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < treshold),]
	downParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] < treshold),]
	upRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(upParental)),]
	downRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(downParental)),]
	cat("--------2---------\n")
	for(x in rownames(upRils)){
		cur <- splitRow.internal(x, upRils, upParental, overlapInd, proportion, margin, ril$parental$groups, 1)
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,x)
		}
	}
	for(x in rownames(downRils)){
		cur <- splitRow.internal(x, downRils, downParental, overlapInd, proportion, margin, ril$parental$groups, 0)
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,x)
		}
	}
	cat("--------3---------\n")
	cat(dim(output),markerNames,"\n")
	ril$rils$genotypes$simulated <- output
	colnames(ril$rils$genotypes$simulated) <- colnames(upRils)
	rownames(ril$rils$genotypes$simulated) <- markerNames
	invisible(ril)
}

splitRow.internal <- function(x, rils, parental, overlapInd, proportion, margin, groupLabels, up=1){
	result <- rep(0,length(rils[x,]))
	A <- parental[which(rownames(parental) == x),which(groupLabels==0)]
	B <- parental[which(rownames(parental) == x),which(groupLabels==1)]
	splitVal <- mean(mean(A,na.rm=TRUE),mean(B,na.rm=TRUE))
	a <- rils[x,] > splitVal
	b <- rils[x,] < splitVal
	if(length(proportion)==2){
		if(up==1){
			genotypes <- c(0,1)
		}else if(up==0){
			genotypes <- c(1,0)
		}
		result[which(a)] <- genotypes[1]
		result[which(b)] <- genotypes[2]
		result[which(rils[x,] == splitVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}else if(length(proportion)==3){
		if(up==1){
			genotypes <- c(0,1,2)
		}else if(up==0){
			genotypes <- c(2,1,0)
		}
		subSplitValDown <- mean(b)
		subSplitValUp <- splitVal + (splitVal - mean(b))
		result[which(rils[x,] < subSplitValDown )] <- genotypes[1]
		result[which(rils[x,] > subSplitValUp )] <- genotypes[3]
		result[which((rils[x,] > subSplitValDown )&&(rils[x,] < subSplitValUp ))] <- genotypes[2]
		result[which(rils[x,] == splitVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}
	invisible(result)
}

filterRow.internal <- function(result, overlapInd, proportion, margin, genotypes){
	if(length(proportion)==2){
		genotypes2 <- genotypes[c(2,1)]
	}else if(length(proportion)==3){
		genotypes2 <- genotypes[c(3,2,1)]
	}
	if(filterRowSub.internal(result,overlapInd,proportion, margin, genotypes)){
		invisible(result)
	}else if(filterRowSub.internal(result,overlapInd,proportion, margin, genotypes)){
		invisible(result)
	}else{
		return(NULL)
	}
}

filterRowSub.internal <- function(result,overlapInd,proportion, margin, genotypes){
	for(i in 1:length(proportion)){
		cur <- sum(result==genotypes[i])/length(result) * 100
		upLimit <- proportion[i] + margin/2
		downLimit <- proportion[i] - margin/2
		#cat(result,"\n",cur,upLimit,downLimit,"\n")
		if(!((cur < upLimit)&&(cur > downLimit))){
			return(FALSE)
		}
	}
	return(TRUE)
}
