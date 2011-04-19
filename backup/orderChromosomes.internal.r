
############################################################################################################
#removeChromosomes.internal: subfunction of filterGenotypes.internal, filtering one row
# 
# cross - object of R/qtl cross type
# minChrLength -if maximal distance between the markers in the chromosome is lower than this value,
#	whole chromosome will be dropped
#
############################################################################################################
orderChromosomes.internal <- function(cross){
	map <- cross$maps$physical[[1]]
	chrtable <- table(map[,1])
	result <- matrix(0,nchr(cross),length(chrtable))
	rownames(result) <- chrnames(cross)
	colnames(result) <- names(chrtable)
	genes <- matrix(1,sum(nmar(cross)),1)
	rownames(genes) <- markernames(cross)
	genes[,1] <- rep(1:length(nmar(cross)),nmar(cross))
	for(i in markernames(cross)){
		chr1 <- genes[i,1]
		chr2 <- map[which(rownames(map)==i),1]
		result[chr1,chr2] <- result[chr1,chr2] + 1
	}
	for(i in 1:length(table(genes[,1]))){
		if(result[i,
	}
}