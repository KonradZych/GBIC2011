# Basic QTL Mapping functions
# 
# (c) 2011-2015 Konrad Zynch
# Version: 0.0.1
#
# Contains:
# qtlbyttest - Basic single marker mapping by using a t.test statistic (For RIL)

#Loading data

qtlbyttest <- function(phenotypes,genotypes,trait){
	res<-NULL
	for(m in 1:ncol(genotypes)){
		pheno_class1 <- phenotypes[which(genotypes[,m]==1),trait]
		pheno_class2 <- phenotypes[which(genotypes[,m]==2),trait]
		res <- c(res,-log10(t.test(pheno_class1,pheno_class2)[[3]]))
	}
	res
}

qtlbyttest_test <- function(){
	setwd("d:/data")
	phenotypes <- read.table()
	genotypes <- read.table()
	res <- qtlbyttest(phenotypes,genotypes,1)
	if(res[20]==0.044){
		cat("test succesfull")
	}else{
		stop("Error wrong answer from qtlbyttest: ",res[20]," should be ",0.044)
	}
}

#qtlbyttest_test()