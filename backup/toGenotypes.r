############################################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified May 2011
# last modified in version: 0.5.1 
# in current version: active, in main workflow
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: toGenotypes
#           convertToGenotypes.internal, splitRow.internal, filterGenotypes.internal, filterRow.internal
#			sortMap.internal, majorityRule.internal, mergeChromosomes.internal, splitRowSubEM.internal
#			checkMu.internal, removeChromosomes.internal
#
############################################################################################################

############################################################################################################
#toGenotypes: Function that chooses from the matrix only appropriate markers with specified rules
# 
# ril - Ril type object, must contain parental phenotypic data.
# use - Which genotypic matrix should be saved to file, real - supported by user and read from file, 
#	simulated - made by toGenotypes, ap - simulated data orderd using gff map
# treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# minChrLength -if maximal distance between the markers in the chromosome is lower than this value,
#	whole chromosome will be dropped
# ... - Parameters passed to formLinkageGroups.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
toGenotypes <- function(ril, use=c("real","simulated","map"), splitMethod=c("EM","mean"),treshold=0.01, overlapInd = 0, proportion = c(50,50), margin = 15, minChrLength = 0, verbose=FALSE, debugMode=0,...){
	#*******CHECKS*******
	s<-proc.time()
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(any(!(is.numeric(ril$parental$phenotypes)))){
		ril <- intoRil(ril, ril$parental$phenotypes)
	}
	if(any(!(is.numeric(ril$rils$phenotypes)))){
		ril <- intoRil(ril, children=ril$rils$phenotypes)
	}
	if(overlapInd < 0 || overlapInd > ncol(ril$rils$phenotypes)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")
	
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	s1 <- proc.time()
	ril <- convertToGenotypes.internal(ril, splitMethod, treshold, overlapInd, proportion, margin, verbose, debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
	
	
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
	cross <- genotypesToCross.internal(ril,use=use,verbose=verbose,debugMode=debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Creating cross object done in:",(e1-s1)[3],"seconds.\n")
	
	#*******ENHANCING CROSS OBJECT*******
	if(use!="map"){
		### FormLinkage groups
		cross <- invisible(formLinkageGroups(cross,reorgMarkers=TRUE,verbose=verbose,...))
		
		### remove shitty chromosomes
		cross <- removeChromosomes.internal(cross,minChrLength)
		### saving as separated object, beacause orderMarkers will remove it from cross object
		removed <- cross$rmv
		
		### Order markers
		#cross <- orderMarkers(cross, use.ripple=TRUE, verbose=verbose)
		
		### Adding real maps		
		if(!(is.null(ril$rils$map))){
			ril <- sortMap.internal(ril)
			cross$maps$physical <- ril$rils$map
			### Majority rule used to order linkage groups
			cross <- segregateChromosomes.internal(cross)
		}
		
		### adding info about removed markers
		cross$rmv <- removed
	}
	
	#*******RETURNING CROSS OBJECT*******
	e<-proc.time()
	if(verbose) cat("toGenotypes done in",(e-s)[3],"seconds\n")
	invisible(cross)
}

############################################################################################################
#convertToGenotypes.internal: function splitting differentially expressed markers into two genotypes
# 
# ril - Ril type object, must contain parental phenotypic data.
# treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information 
#
############################################################################################################
convertToGenotypes.internal <- function(ril, splitMethod, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
	### initialization
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	output <- NULL
	markerNames <- NULL 
	
	### selection step
	upParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < treshold),]
	downParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] < treshold),]
	upRils <- ril$rils$phenotypes[rownames(upParental),]
	downRils <- ril$rils$phenotypes[rownames(downParental),]
	
	### checking if anything is selected and if yes - processing
	if(!(is.null(dim(upRils)))){
		if(!(is.null(dim(downRils)))){
			# best situation
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers and ",nrow(downRils),"downregulated markers.\n")
			cur <- splitRow.internal(downRils, downParental, splitMethod, overlapInd, proportion, margin, ril$parental$groups, 0)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers.\n")
		}
		cur <- splitRow.internal(upRils, upParental, splitMethod, overlapInd, proportion, margin, ril$parental$groups, 1)
		output <- rbind(output,cur[[1]])
		markerNames <- c(markerNames,cur[[2]])
	}else{
		if(!(is.null(dim(downRils)))){
			if(verbose) cat("Selected ",nrow(downRils),"downregulated markers.\n")
			cur <- splitRow.internal(downRils, downParental, splitMethod, overlapInd, proportion, margin, ril$parental$groups, 0)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			stop("None of the markers was selected using specified treshold: ",treshold,"\n")
		}
	}
	
	### putting results inside ril object
	if(is.null(dim(output))) stop("No markers selected.")
	ril$rils$genotypes$simulated <- output
	colnames(ril$rils$genotypes$simulated) <- colnames(upRils)
	rownames(ril$rils$genotypes$simulated) <- markerNames
	invisible(ril)
}

############################################################################################################
#splitRow.internal: subfunction of convertToGenotypes.internal, splitting children markers using parental
# mean values
# 
# rils - matrix of up/down regulated genes in rils
# parental - matrix of up/down regulated genes in parents
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# up - 1 - genes up 0 - down regulated
#
############################################################################################################
splitRow.internal <- function(rils, parental, splitMethod, overlapInd, proportion, margin, groupLabels, up){
	output <- NULL
	markerNames <- NULL
	cat("ril:",nrow(rils),"parental:",nrow(parental),"\n")
	for(x in rownames(rils)){
		if(splitMethod=="mean"){
			cur <- splitRowSub.internal(x, rils, parental, overlapInd, proportion, margin, groupLabels, up)
		}else if(splitMethod=="EM"){
			cur <- splitRowSubEM.internal(x, rils, parental, overlapInd, proportion, margin, groupLabels, up)
		}
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,x)
		}
	}
	invisible(list(output,markerNames))
}

############################################################################################################
#splitRowSub.internal: subfunction of splitRow.internal, splitting one row
# 
# x - name of currently processed row
# rils - matrix of up/down regulated genes in rils
# parental - matrix of up/down regulated genes in parents
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# up - 1 - genes up 0 - down regulated
#
############################################################################################################
splitRowSub.internal <- function(x, rils, parental, overlapInd, proportion, margin, groupLabels, up=1){
	### initialization
	result <- rep(0,length(rils[x,]))
	
	### splitting
	if(length(proportion)==2){
		### splitting into 0/1 genotypes
		if(up==1){
			genotypes <- c(0,1)
		}else if(up==0){
			genotypes <- c(1,0)
		}
		A <- parental[x,which(groupLabels==0)]
		B <- parental[x,which(groupLabels==1)]
		splitVal <- mean(mean(A,na.rm=TRUE),mean(B,na.rm=TRUE))
		rils[x,which(is.na(rils[x,]))] <- splitVal
		a <- rils[x,] > splitVal
		b <- rils[x,] < splitVal
		result[which(a)] <- genotypes[1]
		result[which(b)] <- genotypes[2]
		result[which(rils[x,] == splitVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}else if(length(proportion)==3){
		### splitting into 0/1/2 genotypes
		if(up==1){
			genotypes <- c(0,1,2)
		}else if(up==0){
			genotypes <- c(2,1,0)
		}
		meanVal <- mean(rils[x,],na.rm=TRUE)
		sdVal <- sd(rils[x,],na.rm=TRUE)
		subSplitValDown <- meanVal - sdVal
		subSplitValUp <- meanVal + sdVal
		result[which(rils[x,] < subSplitValDown )] <- genotypes[1]
		result[which(rils[x,] > subSplitValUp )] <- genotypes[3]
		result[which((rils[x,] > subSplitValDown )&&(rils[x,] < subSplitValUp ))] <- genotypes[2]
		result[which(rils[x,] == meanVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}
	
	invisible(result)
}

############################################################################################################
#filterRow.internal : removing from genotypic matrix genes that are no passing specified requirments
# 
# ril - Ril type object, must contain parental phenotypic data.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information

############################################################################################################
filterRow.internal <- function(result, overlapInd, proportion, margin, genotypes){
	### creating inverted genotypes matrix, to be sure, that we won't filter out anythng in correct proportion
	### this function returns either unchanged result vector, which is then rbinded to other results, or
	### NULL, which is ignored by rbind, and we drop current result
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

############################################################################################################
#filterRowSub.internal : subfunction of filterGenotypes.internal, filtering one row
# 
# genotypeRow - currently processed row
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
#
############################################################################################################
filterRowSub.internal <- function(result, overlapInd, proportion, margin, genotypes){
	for(i in 1:length(proportion)){
		cur <- sum(result==genotypes[i])/length(result) * 100
		if(is.na(cur)) return(FALSE)
		upLimit <- proportion[i] + margin/2
		downLimit <- proportion[i] - margin/2
		#cat(result,"\n",cur,upLimit,downLimit,"\n")
		if(!((cur < upLimit)&&(cur > downLimit))){
			return(FALSE)
		}
	}
	return(TRUE)
}

############################################################################################################
#sortMap.internal: subfunction of filterGenotypes.internal, filtering one row
# 
# ril - Ril type object
#
############################################################################################################
sortMap.internal <- function(ril){
	genes <- ril$rils$map
	result <- NULL
	nchr <- length(table(genes[,1]))
	lengths <- vector(mode="numeric",length=nchr+1)
	lengths[1] <- 0
	for(i in 1:nchr){
		current <- genes[which(genes[,1]==i),]
		lengths[i+1] <- max(current[,2]) + lengths[i] + 20000
		current <- current[names(sort(current[,2])),]
		current[,2] <- current[,2] + lengths[i]
		result <- rbind(result,current)
	}
	ril$rils$map <- list(result,lengths[-length(lengths)])
	invisible(ril)
}

############################################################################################################
#segragateChromosomes.internal  - ordering chromosomes using physical map
# 
# cross - object of R/qtl cross type
#
############################################################################################################
segregateChromosomes.internal <- function(cross){
	if(is.null(cross$maps$physical)){
		cat("WARNING: no physical map, function will return unchanged cross object\n")
	}else{
		output <- majorityRule.internal(cross)
		print(output)
		### until every chr on phys map is match exactly once
		while(max(apply(output,2,sum))>1){
			toMerge <- which(apply(output,2,sum)>1)
			for(curToMerge in toMerge){
				curMerge <- which(output[,curToMerge]==max(output[,curToMerge]))
				map <- cross$maps$physical
				cross <- mergeChromosomes.internal(cross,curMerge,curMerge[1])
				cross$maps$physical <- map
				output <- majorityRule.internal(cross)
			}
		}
		
		order1 <- matrix(0,ncol(output),nrow(output))
		order2 <- matrix(1,ncol(output),nrow(output))
		### until next iteration doesn't change the result
		while(any(order1!=order2)){
			order1 <- output
			for(l in 1:ncol(output)){
				cur <- which(output[,l]==max(output[,l]))
				if(cur!=l)cross <- switchChromosomes.internal(cross,cur,l)
				output <- majorityRule.internal(cross)
			}
			order2 <- output
		}
		names(cross$geno) <- 1:length(cross$geno)
	}
	invisible(cross)
}

############################################################################################################
#majorityRule.internal - subfunction of segragateChromosomes.internal, returns matrix showing for every
# reco map chromosome from which physicall map chromosome majority of markers comes
# 
# cross - object of R/qtl cross type
#
############################################################################################################
majorityRule.internal <- function(cross){
	knchrom <- length(table(cross$maps$physical[[1]][,1]))
	result <- matrix(0, length(cross$geno), knchrom)
	output <- matrix(0, length(cross$geno), knchrom)
	for(i in 1:length(cross$geno)){
		cur_ys <- colnames(cross$geno[[i]]$data)
		cur_xs <- cross$maps$physical[[1]][cur_ys,]
		for(j in 1:knchrom){
			result[i,j] <- sum(cur_xs[,1]==j)/nrow(cur_xs)
		}
		output[i,which(result[i,]==max(result[i,]))] <- 1
	}
	rownames(result) <- 1:nrow(result)
	colnames(result) <- 1:ncol(result)
	rownames(output) <- 1:nrow(output)
	colnames(output) <- 1:ncol(output)
	
	if(min(apply(output,2,max))!=1){
		toCheck <- which(apply(output,2,sum)!=1)
		for(x in toCheck){
			output[,x] <- 0
			output[which(result[,x]==max(result[,x])),x] <- 1
		}
	}	
	invisible(output)
}

############################################################################################################
#mergeChromosomes.internal - subfunction of segragateChromosomes.internal, merging multiple chromosomes into
# one
# 
# cross - object of R/qtl cross type
# chromosomes - chromosomes to be merged
# name - chromosome
#
############################################################################################################
mergeChromosomes.internal <- function(cross, chromosomes, name){
	cat("Merging chromosomes",chromosomes,"to form chromosome",name,"names:",names(cross$geno),"\n")
	geno <- cross$geno
	markerNames <- NULL
	for(j in chromosomes){
		if(j!=name) markerNames <- c(markerNames, colnames(geno[[j]]$data))
	}
	for(k in markerNames) cross <- movemarker(cross, k, name)
	cat("Ordering markers on newly merged chromosome\n")
	#cross <- orderMarkers(cross, chr=name)
	invisible(cross)
}

############################################################################################################
#splitRowSubEM.internal - subfunction of segragateChromosomes.internal, merging multiple chromosomes into
# one
# 
# cross - object of R/qtl cross type
# chromosomes - chromosomes to be merged
# name - chromosome
#
############################################################################################################
splitRowSubEM.internal <- function(x, rils, parental, overlapInd, proportion, margin, groupLabels, up=1){
	### initialization
	print(x)
	nrDistributions <- length(proportion)
	result <- rep(0,length(rils[x,]))
	EM <- cat(normalmixEM(sort(rils[x,]), k=nrDistributions, maxrestarts=0, maxit = 100,fast=TRUE),file=NULL)
	if(up==1){
		genotypes <- c(0:(nrDistributions-1))
	}else if(up==0){
		genotypes <- c((nrDistributions-1):0)
	}
	len <- vector(mode="numeric",length=nrDistributions)
	for(i in 1:nrDistributions){
		len[i]<-length(rils[x,])*EM$lambda[i]
		startVal <- sum(len[1:i-1])
		rils[x,which(rils[x,] %in% sort(rils[x,])[startVal:(startVal+len[i])])] <- genotypes[i]
	}
	#result <- checkMu.internal(EM)
	if(checkMu.internal(EM)){
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}else{
		result<- NULL
	}
	invisible(result)
}

############################################################################################################
#checkMu.internal: checking if fitted normal distributions do not overlap
# 
# EM - output of normalmixEM function
#
############################################################################################################
checkMu.internal <- function(EM){
	for(i in 2:length(EM$mu)){
		if((EM$mu[i]-EM$sigma[i])<(EM$mu[i-1]+EM$sigma[i-1])) return(FALSE)
	}
	return(TRUE)
}
