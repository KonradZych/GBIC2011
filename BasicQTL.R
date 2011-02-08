# Basic QTL Mapping functions
# 
# (c) 2011-2015 Konrad Zych
# Version: 0.0.1
#
# Contains:
#qtlbyttest - Basic single marker mapping by using a t.test statistic (For RIL)
#heatmapqtl - Creates data file concerning all traits and also image of it

#qtlbyttest
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#trait: integer value of the column to analyse
#return: Vector of significance (-LOD10) of linkage to that marker
qtlbyttest <- function(phenotypes,genotypes,trait){
	res<-NULL
	for(m in 1:ncol(genotypes)){
		pheno_class1 <- phenotypes[which(genotypes[,m]==1),trait]
		pheno_class2 <- phenotypes[which(genotypes[,m]==2),trait]
		res <- c(res,-log10(t.test(pheno_class1,pheno_class2)[[3]]))
	}
	res
}

#heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#result: name of the resulting array
#return: translated matrix of results of qtlbyttest and an image of it
heatmapqtl <- function(phenotypes,genotypes,output){
	output <- NULL
	for(y in 1:ncol(phenotypes)){
		output <-rbind(output,qtlmap(phenotypes,genotypes,y))
	}
	output <- t(output)
	image(output)
	output
}


qtlbyttest_test <- function(){
	#Creating vector of correct data - O RLY?!
	cor<-c(0.32874301, 0.54131594, 0.78309899, 1.38776104, 0.52714749,
	0.41091438, 4.46812544, 1.19857837, 1.51278035, 1.23372805, 0.06181733,
	0.24965065, 1.65847530, 2.73929366, 0.68674918, 3.84755060, 1.89903764,
	0.77890058, 1.03957367, 0.30353727, 0.11931368, 0.83799584, 0.10603482, 0.43340857)
	setwd("d:/data")
	#Loading data
	phenotypes <- as.matrix(read.table("measurements_ordered.txt", sep=""))
	genotypes <- as.matrix(read.table("genotypes_ordered.txt", sep=""))
	#Testing loop
	for(x in 1:ncol(genotypes)){
		res <- qtlbyttest(phenotypes,genotypes,x)
		if(res[116]==cor[x]){
			cat("Value for the trait ",x,"is correct.") 
		}else{
			stop("Error wrong value for trait ",x," is ",res[116]," should be ",cor[x])
		}
		cat("All values are correct.")
	}
}

#qtlbyttest_test()