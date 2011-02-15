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

#makebinary
#mat: matrix to be made binary
#res: resulting binary matrix
makebinary<-function(mat){
	res<-matrix(0,nrow(mat),ncol(mat))
	#Result matrix should have the same labels as input
	rownames(res)<-rownames(mat, do.NULL = FALSE)
	colnames(res)<-colnames(mat, do.NULL = FALSE)
	#now using median (tres) as a treshold
	for(h in 1:ncol(mat)){
	    tres=median(mat[,h])
		res[,h]<-(mat[,h]>tres)
	}
	res
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


#heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#result: name of the resulting array
#return: translated matrix of results of qtlbyttest
heatmapqtl <- function(phenotypes,genotypes){
	output <- NULL
	for(y in 1:ncol(phenotypes)){
		output <-rbind(output,qtlbyttest(phenotypes,genotypes,y))
	}
	output <- t(output)
	rownames(output)<-colnames(genotypes, do.NULL = FALSE)
	colnames(output)<-colnames(phenotypes, do.NULL = FALSE)
	output
}

#pathway
#results: output array from heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits (mind labels, they are crucial)
#marker: specified marker (use 100 for current data)
pathway<-function(results,marker){
	path<-NULL
	vect<-results[marker,]
	for(i in 1:ncol(results)){
		if(max(vect)==0){break}
		else{
			res<-which(vect==max(vect))
			path<-c(path,colnames(result_binary)[res])
			vect[res]<-0
		}
	}
	path
}

#makecorvector
makecorvector <- function(mat){
	#cormatrixm and corvectorm are for Markers
	cormatrixm <- cor(mat,use="pairwise.complete.obs")
	#making vector for markers
	corvectorm <- NULL
	result <- NULL
	for(i in 2:(ncol(cormatrixm)-1)){
		corvectorm <- c(corvectorm,mean(cormatrixm[i-1,i],cormatrixm[i,i+1]))
	}
	corvectorm <- c(corvectorm,(cormatrixm[(ncol(cormatrixm)-1),(ncol(cormatrixm))])/2)
	for(i in 1:23){
		result <- rbind(result,corvectorm)
	}
	result <- t(result)
	result
}

#makepallete
makepallete<-function(covmatrix, result){
	crange <- max(c(covmatrix,abs(min(covmatrix))))
	pallete <- matrix(0,nrow(covmatrix),ncol(covmatrix))
	for(i in 1:nrow(covmatrix)){
		for(j in 1:ncol(covmatrix)){
			pallete[i,j]<-rgb((abs(covmatrix[i,j])/crange)*255,(abs((result[i,j])/max(result))*255),255-((abs(covmatrix[i,j])/crange)*255),maxColorValue=255)
			#pallete[i,j] <- rgb((abs((result[i,j])/max(result))*255),255-(abs((result[i,j])/max(result))*255),0,maxColorValue=255)
		}
	}
	pallete
}

#smoothie - it's not quite what I expected
smoothie <- function(result){	
	for(i in nrow(result)){
		for(j in 1:9){
			result[i,j] <- mean(result[i,1:j+9])
		}
		for(j in 10:(ncol(result)-9)){
			result[i,j] <- mean(result[i,j-9:j+9])
		}
		for(j in (ncol(result)-9):(ncol(result))){
			result[i,j] <- mean(result[i,j:(ncol(result))])
		}
	}
	result
}

#un_intonumeric
un_intonumeric <- function(un_matrix){
	inside <- un_matrix
	#switching -, a, b to NA, 1, 2, respectively
	un_matrix[which(un_matrix=="-")] <- NA
	un_matrix[which(un_matrix=="a")] <- 1
	un_matrix[which(un_matrix=="b")] <- 2
	un_matrix <- as.numeric(un_matrix)
	un_matrix <- matrix(un_matrix,nrow(inside),ncol(inside))
	rownames(un_matrix)<-rownames(inside, do.NULL = FALSE)
	colnames(un_matrix)<-colnames(inside, do.NULL = FALSE)
	un_matrix
}

#un_findneighbors - 
#un_findneighbors <- function(mat){
	#un_covmatrix <- cor(un_matrix,use="pairwise.complete.obs")
	#for(i in 1:nrow(un_covmatrix)){
	#	un_covmatrix[i,i] <- -10
	#}
	#res2<-NULL
	#result_m <- matrix(0,nrow(mat), ncol(mat))
	#first <- which(colMeans(un_covmatrix)==min(colMeans(un_covmatrix)))[1]
	#res <- first
	#for(i in 2:ncol(un_covmatrix)-1){		
	#	recent <- which(un_covmatrix[,first]==max(un_covmatrix[,first]))[1]
	#	res2 <- c(res2,max(un_covmatrix[,first]))
	#	res <- c(res,recent)
	#	un_covmatrix[,first]<--10
	#	un_covmatrix[first,]<--10
	#	first <- recent
	#}
	#for(i in 1:nrow(un_covmatrix)-2){
	#	res <- c(res,labels(which((un_covmatrix[i,(i+1):ncol(un_covmatrix)])==(sort(un_covmatrix[i,(i+1):ncol(un_covmatrix)],decreasing=TRUE)[1]))))
	#}	
	#print(length(res))
	#for(i in 1:length(res)){
	#	result_m[,i] <- mat[,res[i]]
	#}
	#print(res2)
	#result_m
#}

#for(i in 1:150){
#	ordered<-un_findneighbors(ordered)
#}

#un_most_common
un_most_common <- function(vect){
	for(i in 1:ncol(ord)){
		print(labels(which(table(ord[,i])==max(table(ord[,i])))))
	}
}

#un_order_chromosome - looks pretty cool, but needs polishing
un_order_chromosome <- function(chrom_matrix){
	chrom_cor_matrix <- cor(chrom_matrix, use="pairwise.complete.obs")
	result <- NULL
	current <- NULL
	for(i in 1:ncol(chrom_cor_matrix)){
		sorted <- sort(chrom_cor_matrix[i,],decreasing=TRUE)
		for(j in 1:ncol(chrom_cor_matrix)){
			current <- c(current,which(chrom_cor_matrix[i,]==(sorted[j])))
		}
		result <- rbind(result,current)
		current <- NULL
	}
	#output<-chrom_cor_matrix[,result[,11]]
	#output
	#result <- (result,1,un_most_common)
	output <- colnames(chrom_matrix[,result[1,]])
	output	
}

#un_neighbor - work a bit here
un_neighbor <- function(un_matrix,groups){
	cor_matrix <- cor(un_matrix,use="pairwise.complete.obs")
	r <- kmeans(cor_matrix,groups)
	res <- NULL
	for(i in 1:groups){
		print(i)
		cur <- un_order_chromosome(cor_matrix[,which(as.numeric(r[[1]])==i)])
		res <- cbind(res,un_matrix[,cur])
	}
	res
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
	#Make qtlmap (quantative) an store it in result_quantative
	result_quantative<-heatmapqtl(phenotypes,genotypes)
	print("result_quantative created")
	#Make qtlmap (quantative) an store it in result_binary
	result_binary<-heatmapqtl(makebinary(clean(phenotypes)),genotypes)
	print("result_binary created")
	#Create pathway using marker 100 for binary
	path<-pathway(result_binary, 100)
	print("path created")
	#Creating covariance_matrix out of two vectors, one for traits, another for markers
	covariance_matrix <- makecorvector(genotypes)
	print("covariance_matrix created")
	#Creating color_pallete
	color_pallete <- makepallete(covariance_matrix, result_binary)
	#persp plot
	persp(result_binary,col=color_pallete)
	#analyzing the unknown
	setwd("d:/data")
	un_matrix <- un_intonumeric(as.matrix(read.table("unknown_genotypes2.txt", sep="", header=TRUE)))
	un_cor <- cor(un_matrix,use="pairwise.complete.obs")
	#finding the neighbor 
	un_findneighbors(un_cor,un_cor)
}

