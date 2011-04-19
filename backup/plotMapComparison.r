	nchrom1 <- length(cross1$geno)
	chrom1 <- chrnames(cross1)
	nchromY <- length(cross2$geno)
	chrom2 <- chrnames(cross2)
	nchrom <- (chrom1 %in% chrom2)
	cat("Object cross1 contains chromosomes:",paste(chrom1,sep=", "),"and cross2 contains chromosomes:",paste(chrom2,sep=", "),".\n")
	nchromX <- 1:nchrom1
	if(!is.null(chr)){
		nchromX <- NULL
		while(length(chr)>=1){
			if(!(chr[1] %in% nchrom)){
				warning("There are only chromosomes: ",paste(nchrom,sep=", ")," chromosome: ",chr[1]," not found.\n")
			}else{
				nchromX <- c(nchromX,chr[1])
			}
			chr <- chr[-1]
		}
	}
	plotChromosomeMap.internal(cross1,cross2,nchromX,nchromY)

}

compareGeneLocation.internal <- function(cross1, cross2){
	genes1 <- compareGeneLocationSub.internal(cross1)
	genes2 <- compareGeneLocationSub.internal(cross2)
	genes1 <- mapMarkers.internal(genes1,genes2,mapMode=1)
	genes2 <- mapMarkers.internal(genes2,genes1,mapMode=1)
	result <- matrix(0,nrow(genes1),6)
	rownames(result) <- rownames(genes1)
	for(i in rownames(result)){
		result[i,c(1,2,3)] <- genes1[i,c(1,2,3)]
		result[i,c(4,5,6)] <- genes2[i,c(1,2,3)]
	}
	invisible(result)
}

compareGeneLocationSub.internal <- function(cross){
	genes <- matrix(1,sum(nmar(cross)),3)
	rownames(genes) <- markernames(cross)
	genes[,1] <- rep(1:length(nmar(cross)),nmar(cross))
	genes <- genePosition.internal(genes,cross)
	invisible(genes)
}

genePosition.internal <- function(genes,cross){
	for(i in rownames(genes)){
		chr <- genes[i,1]
		if(chr>1){ 
			prev <- sum(chrlen(cross)[1:(chr-1)]) + (0.13* max(chrlen(cross)) * (chr-1))
		}else{ 
			prev <- 0
		}
		genes[i,2] <- pull.map(cross)[[chr]][i]
		genes[i,3] <- genes[i,2] + prev
	}
	invisible(genes)
}

plotChromosomeMap.internal <- function(cross1,cross2,nchromX,nchromY){
	genes <- compareGeneLocation.internal(cross1,cross2)
	genes <- genes[which(genes[,1] %in% nchromX),]
	plot(x=genes[1,3], y=genes[1,6], xlim=c(min(genes[,3]),max(genes[,3])), ylim=c(min(genes[,6]),max(genes[,6])), col="red", xlab="Reference cross", ylab="Cross 2", main="Comparison of genetic maps")
	if(max(nchromY) < max(nchromX)){
		cl <- topo.colors(max(nchromX))
	}else{
		cl <- c(topo.colors(max(nchromX)),rep("red",(max(nchromY) - max(nchromX))))
	}
	for(i in rownames(genes)){
		mark <- 0 + genes[i,1]
		color <- cl[genes[i,4]]
		points(x=genes[i,3], y=genes[i,6],col=color,pch=mark)
	}
	for(x in nchromX){
		if(x > 1) abline(v=sum(chrlen(cross1)[1:(x-1)]) + (0.13* max(chrlen(cross1)) * (x-1)),lty=2)
	}
	for(y in 1:nchromY){
		if(y > 1) abline(h=sum(chrlen(cross2)[1:(y-1)]) + (0.13* max(chrlen(cross2)) * (y-1)),lty=2)
	}
}