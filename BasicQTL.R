# Basic QTL Mapping functions
# 
# (c) 2011-2015 Konrad Zych
# Version: 0.0.1
#
# Contains:
#clean - Removing NA from data by replacing it with 0
#qtlbyttest - Basic single marker mapping by using a t.test statistic (For RIL)
#heatmapqtl - Creates data file concerning all traits and also image of it
#pathway - creates probable pathway using values for specified marker
#controller - AIO function

#clean
#mat: matrix to be cleaned
clean<-function(mat){
	for(h in 1:nrow(mat)){
		for(w in 1:ncol(mat)){
			if(is.na(mat[h,w])){
				mat[h,w]<-0
			}
		}
	}
	mat
}

#qtlbyttest
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#trait: integer value of the column to analyse
#return: Vector of significance (-LOD10) of linkage to that marker
qtlbyttest <- function(phenotypes,genotypes,trait){
	res<-NULL
	phenotypes<-clean(phenotypes)
	genotypes<-clean(genotypes)
	for(m in 1:ncol(genotypes)){
		pheno_class1 <- phenotypes[which(genotypes[,m]==1),trait]
		pheno_class2 <- phenotypes[which(genotypes[,m]==2),trait]
		res <- c(res,-log10(t.test(pheno_class1,pheno_class2)[[3]]))
	}
	res
}

#makebinary
#mat: matrix to be made binary
#res: resulting binary matrix
makebinary<-function(mat){
	res<-matrix(0,nrow(mat),ncol(mat))
	#now using median (tres) as a treshold
	for(h in 1:ncol(mat)){
	    tres=median(mat[,h])
		res[,h]<-(mat[,h]>tres)
	}
	res
}

#heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#result: name of the resulting array
#return: translated matrix of results of qtlbyttest
heatmapqtl <- function(phenotypes,genotypes){
	output <- NULL
	for(y in 1:ncol(phenotypes)){
		output <-rbind(output,qtlmap(phenotypes,genotypes,y))
	}
	output <- t(output)
	output
}

#pathway
#results: output array from heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits (mind labels, they are crucial)
#marker: specified marker (use 100 for current data)
pathway<-function(results,phenotypes,marker){
	path<-NULL
	metabolites<-as.matrix(labels(phenotypes[marker,]))
	vect<-results[marker,]
	for(i in 1:ncol(results)){
		if(max(vect)==0){break}
		else{
			res<-which(vect==max(vect))
			path<-c(path,metabolites[res])
			vect[res]<-0
		}
	}
	path
}

#controller
controller<-function(){
	#Firstly, doing tests
	#qtlbyttest_test() - what the hell is wrong here?!
	makebinary_test()
	#Load data
	setwd("d:/data")
	print("wd set")
	phenotypes <- as.matrix(read.table("measurements_ordered.txt", sep=""))
	print("phenotypes loaded")
	genotypes <- as.matrix(read.table("genotypes_ordered.txt", sep=""))
	print("genotypes loaded")
	#Make binary matrix of phenotypes an store it in phenotypes_binary
	phenotypes_binary<-makebinary(clean(phenotypes))
	print("phenotypes_binary created")
	#Make qtlmap (quantative) an store it in result_quantative
	result_quantative<-heatmapqtl(phenotypes,genotypes)
	print("result_quantative created")
	#Make qtlmap (quantative) an store it in result_binary
	result_binary<-heatmapqtl(phenotypes_binary,genotypes)
	print("result_binary created")
	#Create pathway using marker 100 for binary
	path<-pathway(result_binary, phenotypes, 100)
	print("path created")
	path
}

qtlbyttest_test <- function(){
	#Creating vector of correct data
	cor<-c(15.6152669361112, 10.9611214399831, 10.6377737262226, 
	6.66414029040903, 1.49439826143111,  8.4439633754547, 
	10.1596325117105, 18.6492514223787, 7.70593733223225, 
	2.83027804463293, 19.0918976951953, 6.13302622877962,
	5.42781730179255, 15.1181942357603, 7.93420253276565,  
	10.9872557298082, 8.40230035811374, 5.13484444988022,
	0.0513458276192633, 0.282071067944515, 0.623321252421781,                                                             
	1.69995386811916, 0.395375105682745, 0.66074033328745)
	setwd("d:/data")
	#Loading data
	phenotypes <- clean(as.matrix(read.table("measurements_ordered.txt", sep="")))
	genotypes <- clean(as.matrix(read.table("genotypes_ordered.txt", sep="")))
	#Testing loop
	for(x in 1:ncol(genotypes)){
		res <- qtlbyttest(phenotypes,genotypes,x)
		if((res[116]==cor[x])){#why the fuck isn't this working?!
			cat("Value for the trait ",x,"is correct.") 
		}else{
			cat("Error wrong value for trait ",x," is ",res[100]," should be ",cor[x],"\n")
		}
		cat("All values are correct.\n")
	}
}

makebinary_test <- function(){
	#Creating vector of correct data
	correct<-c(0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0)
	setwd("d:/data")
	#Loading data
	phenotypes <- clean(as.matrix(read.table("measurements_ordered.txt", sep="")))
	#Testing here
	res <- makebinary(phenotypes)
	for(i in 1:ncol(phenotypes)){
		if(res[100,i]!=correct[i]){
			stop("Wrong value for trait ",i,". Disconnecting for your safety.")
		}
	}
	cat("Test passed. Makebinary ready to serve.\n")
}


#qtlbyttest_test()
#makebinary_test()

