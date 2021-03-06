batcheffectcheck <- function(pheno, cutoff=2, minimumvariance = 0.75){
  correlation <- cor(pheno)
  tree <- heatmap(correlation, Colv=NA, scale="none", keep.dendro=T)
  branches <- cut(tree$Rowv,cutoff)
  returnlist <- vector("list", length(branches$lower))
  for(x in 1:length(branches$lower)){
    returnlist[[x]] <- unlist(branches$lower[[x]])
    names(returnlist[[x]]) <- colnames(pheno)[unlist(branches$lower[[x]])]
    cat("Group",x,":",paste(colnames(pheno)[unlist(branches$lower[[x]])],";"),"\n")
  }
  invisible(returnlist)
}



batcheffectcorrect <- function(pheno, batchlist, minimumvariance = 0.75){
  s <- proc.time()
  cormatrix <- NULL
  cat("",file="tmpbatch.out")
  for(y in 1:ncol(pheno)){
    if(y %% 1000 == 0){
      e <- proc.time()
      cat("Done:",y,"/",ncol(pheno),"in",as.numeric(e[3]-s[3]),"secs\n")
      s <- e
    }
    oamean <- mean(pheno[,y])
    groupmeans <- lapply(lapply(batchlist, fun <- function(x){pheno[x,y]}),mean,na.rm=TRUE)
    diffmeans <- unlist(groupmeans) - oamean
    traitcorrection <- rep(0,nrow(pheno))
    for(x in 1:length(diffmeans)){
      traitcorrection[batchlist[[x]]] <- diffmeans[x]
    }
    cat(traitcorrection,"\n",sep="\t",file="tmpbatch.out",append=TRUE)
  }
  cormatrix <- t(read.table("tmpbatch.out",sep="\t"))
  cormatrix <- cormatrix[1:nrow(pheno),1:ncol(pheno)]
}

rilcorrect <- function(ril,cutoff){
	befc <- batcheffectcheck(tom_c, 0.75, 0)
	befc
	corm <- batcheffectcorrect(t(tom_c), befc, 0)
	invisible(corm)
}