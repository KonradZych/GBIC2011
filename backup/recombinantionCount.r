##################################################################################################
#
# recombinationCount.R
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
# Contains: recombinationCount
#           recombinationCountRow, recombinationCountRowSub, recombinationCountRowFlipSub
#
##################################################################################################

#recombinationCount - counting recobinations needen to go from one matrix to another
#genotypicMatrix - rows: markers, cols: individuals
#flip - specifies whether one of the rows that are being compared should be flipped(1) or not(0)
#verbose - standard
#debugMode - standard
recombinationCount <- function(genotypicMatrix,flip=0,verbose=FALSE,debugMode=0){
	genotypicMatrix <- switchMatrixValues(genotypicMatrix,before=c("A","B"),after=c(0,1))
	res <- apply(genotypicMatrix,1,recombinationCountRow,genotypicMatrix,flip)
	res
}

#recombinationCountRow - using recombinationCountRowSub to apply comparison between every two rows
#genotypicMatrixRow - data for single marker
#flip - specifies whether one of the rows that are being compared should be flipped(1) or not(0)
#verbose - standard
#debugMode - standard
recombinationCountRow <- function(genotypicMatrixRow,genotypicMatrix,flip=0,verbose=FALSE,debugMode=0){
	if(flip==0){output <- apply(genotypicMatrix,1,recombinationCountRowSub,genotypicMatrixRow)}
	if(flip==1){output <- apply(genotypicMatrix,1,recombinationCountRowFlipSub,genotypicMatrixRow)}
	output
}

#recombinationCountRowSub - comparing two rows
#genotypicMatrixRow1 & genotypicMatrixRow2 - data for single marker
recombinationCountRowSub <- function(genotypicMatrixRow1,genotypicMatrixRow2){
	sum((genotypicMatrixRow1)!=(genotypicMatrixRow2))
}

#recombinationCountRowFlipSub - comparing two rows
#genotypicMatrixRow1 & genotypicMatrixRow2 - data for single marker
recombinationCountRowFlipSub <- function(genotypicMatrixRow1,genotypicMatrixRow2){
	sum((genotypicMatrixRow1)!=(1-(genotypicMatrixRow2)))
}

recombinationCount.test(){
	
}