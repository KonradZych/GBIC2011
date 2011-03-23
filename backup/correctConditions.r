#####################################################################
#
# bestClustering.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified March 2011
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
# Contains: bestClustering
#				bestClusteringSub, bestClusteringSubSub
#
#####################################################################

correctConditions <- function(expressionMatrix, groups, iterations=100, verbose=FALSE, debugMode=0){
	#CHECKS
	s <- proc.time()
	if(verbose && debugMode==1) cat("bestClustering starting.\n")
	pointsMatrix <- bestClustering(expressionMatrix, groups, iterations, use="c", verbose=verbose, debugMode=debugMode)
	r <- kmeans(pointsMatrix,groups)
	#correctedMatrix <- apply(expressionMatrix,2,correctConditionsSub,r)
	referenceMean <- mean(expressionMatrix[which(r[[1]]==1)])
	sc <- proc.time()
	for(i in 2:groups){
		groupMean <- mean(expressionMatrix[which(r[[1]]==i)])
		corValue <- groupMean-referenceMean
		print(corValue)
		expressionMatrix[which(r[[1]]==i)] <- expressionMatrix[which(r[[1]]==i)]-corValue
	}
	ec <- proc.time()
	if(verbose && debugMode==2)cat("Normalizing inside groups done in:",(ec-sc)[3],"seconds.\n")
	e <- proc.time()
	if(verbose) cat("correctConditions done in",(e-s)[3],"seconds.\n")
	invisible(expressionMatrix)
}
