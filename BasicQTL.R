# Basic QTL Mapping functions
# 
#(c) 2011-2015 Konrad Zych
# Version: 0.0.1
#
# Contains:
#clean - Removing NA from data matrix by replacing it with 0
#qtlbyttest - Basic single marker mapping by using a t.test statistic (For RIL)
#heatmapqtl - Creates data file concerning all traits and also image of it
#pathway - creates probable pathway using values for specified marker
#controller - AIO function

#clean
#matrix_to_be_cleaned: matrix containing data mixed with NA
#returns: matrix_to_be_cleaned - the same matrix with NAs replaced with 0
clean<-function(matrix_to_be_cleaned){
	for(h in 1:nrow(matrix_to_be_cleaned)){
		for(w in 1:ncol(matrix_to_be_cleaned)){
			if(is.na(matrix_to_be_cleaned[h,w])){
				matrix_to_be_cleaned[h,w]<-0
			}
		}
	}
	matrix_to_be_cleaned
}

#makebinary
#matrix_to_be_made_binary: matrix containing data
#returns: output - matrix containing 0s (value below the treshold) and 1 (value above the treshold) currently, the treshold is median
makebinary<-function(matrix_to_be_made_binary){
	output<-matrix(0,nrow(matrix_to_be_made_binary),ncol(matrix_to_be_made_binary))
	#resulting matrix should have the same labels as input
	rownames(output)<-rownames(matrix_to_be_made_binary, do.NULL = FALSE)
	colnames(output)<-colnames(matrix_to_be_made_binary, do.NULL = FALSE)
	#using median as a treshold value
	for(h in 1:ncol(matrix_to_be_made_binary)){
	    tres=median(matrix_to_be_made_binary[,h])
		output[,h]<-(matrix_to_be_made_binary[,h]>tres)
	}
	output
}

#qtlbyttest
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#trait: integer value of the column to analyse
#return: output - Vector of significance (-LOD10) of linkage to that marker
qtlbyttest <- function(phenotypes,genotypes,trait){
	output<-NULL
	phenotypes<-clean(phenotypes)
	genotypes<-clean(genotypes)
	for(m in 1:ncol(genotypes)){
		pheno_class1 <- phenotypes[which(genotypes[,m]==1),trait]
		pheno_class2 <- phenotypes[which(genotypes[,m]==2),trait]
		output <- c(output,-log10(t.test(pheno_class1,pheno_class2)[[3]]))
	}
	output
}


#heatmapqtl
#phenotypes: Matrix of row: individuals, columns: traits
#genotypes: Matrix of row: individuals, columns: markers
#return: output - translated matrix of results of qtlbyttest
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
#return: output - proposed pathway (linear!) based on height of peak from marker 100
pathway<-function(results,marker){
	output<-NULL
	vect<-results[marker,]
	for(i in 1:ncol(results)){
		if(max(vect)==0){break}
		else{
			output<-which(vect==max(vect))
			output<-c(output,colnames(result_binary)[output])
			vect[output]<-0
		}
	}
	output
}

#makecorvector
#first_matrix: 
#two_matrices:
#second_matrix:
#return: output - 
makecorvector <- function(first_matrix,two_matrices=0,second_matrix=NULL){
	#first_matrix_cor and first_matrix_corv are for Markers
	first_matrix_cor <- cor(first_matrix,use="pairwise.complete.obs")
	#making vector for markers
	first_matrix_corv <- NULL
	output <- NULL
	for(i in 2:(ncol(first_matrix_cor)-1)){
		first_matrix_corv <- c(first_matrix_corv,mean(first_matrix_cor[i-1,i],first_matrix_cor[i,i+1]))
	}
	first_matrix_corv <- c(first_matrix_corv,(first_matrix_cor[(ncol(first_matrix_cor)-1),(ncol(first_matrix_cor))])/2)
	for(i in 1:23){
		output <- rbind(output,first_matrix_corv)
	}
	output <- t(output)
	if(two_matrices==1){
		for(i in 1:nrow(output)){
			for(j in 1:ncol(output)){
				output[i,j]<-output[i,j]+second_matrix[i,j]
			}
		}
	}
	output
}

#make_topo_pallete
make_topo_pallete <- function(result_matrix){
	cur_mean <- mean(result_matrix)
	cur_sd <- mean(sd(result_matrix))
	cur_range <- abs(max(result_matrix)-min(result_matrix))
	topo_pallete<-matrix(0,nrow(result_matrix),ncol(result_matrix))
	for(i in 1:nrow(result_matrix)){
		for(j in 1:ncol(result_matrix)){
			if(result_matrix[i,j]<(cur_mean-cur_sd)){
				topo_pallete[i,j]<-rgb(0,0,abs(result_matrix[i,j]/cur_range)*255,maxColorValue=255)
			}else if((result_matrix[i,j]>(cur_mean+cur_sd))){
				topo_pallete[i,j]<-rgb(abs(result_matrix[i,j]/cur_range)*255,0,0,maxColorValue=255)
			}else{
				topo_pallete[i,j]<-rgb(255-abs(result_matrix[i,j]/cur_range)*255,0,abs(result_matrix[i,j]/cur_range)*255,maxColorValue=255)
			}
			
		}
	}
	topo_pallete
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

#un_best_clustering
un_best_clustering <- function(cor_matrix,nr_iterations){
	res <- NULL
	map <- matrix(0,nrow(cor_matrix),ncol(cor_matrix))
	output <- NULL
	for(i in 1:nr_iterations){
		r <- kmeans(cor_matrix,1)
		res <- rbind(res,(as.numeric(r[[1]])))
	}
	print("OK nr_iterations")
	print(dim(res))
	for(i in 1:nr_iterations){
		for(j in 1){
			for(k in which(res[i,]==j)){
				for(l in which(res[i,]==j)){
					map[k,l] <- map[k,l] + 1
				}
			}
		}
	}
	print("OK map")
	#output<-kmeans(map,10)
	#output
	map
}

#un_order_chromosome 
un_order_chromosome <- function(chrom_matrix){
	chrom_cor_matrix <- cor(chrom_matrix, use="pairwise.complete.obs")
	result <- 1
	current <- NULL
	chrom_cor_matrix[,1]<--10
	for(i in 1:ncol(chrom_cor_matrix)){
		chrom_cor_matrix[i,i]<--10
	}
	i<-1
	while(length(result)<ncol(chrom_cor_matrix)){
		j <- which(chrom_cor_matrix[i,]==max(chrom_cor_matrix[i,]))[1]
		result<-c(result,j)
		chrom_cor_matrix[,j]<--10
		chrom_cor_matrix[j,i]<--10
		i<-j
	}
	output <- colnames(chrom_matrix[,result])
	output	
}

#un_neighbor - work a bit here
un_neighbor <- function(un_matrix,nr_iterations=1000,groups=10){
	cor_matrix <- cor(un_matrix,use="pairwise.complete.obs")
	print("OK cor_matrix")
	r <- un_best_clustering(cor_matrix,nr_iterations)
	print("OK un_best_clustering")
	res <- NULL
	for(i in 1:groups){
		print(i)
		cur <- un_order_chromosome(cor_matrix[,which(r[[1]]==i)])
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
	setwd("D:/data")
	print("wd set")
	phenotypes <- as.matrix(read.table("phenotypes.txt", sep=""))
	print("phenotypes loaded")
	genotypes <- as.matrix(read.table("genotypes.txt", sep=""))
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
	covariance_matrix <- makecorvector(genotypes,1,result_binary)
	print("covariance_matrix created")
	#Creating color_pallete
	color_pallete <- make_topo_pallete(covariance_matrix)
	#persp plot
	persp(result_binary,col=color_pallete)
	#analyzing the unknown
	setwd("D:/data")
	un_matrix <- un_intonumeric(as.matrix(read.table("unknown_genotypes.txt", sep="", header=TRUE)))
	#finding the neighbor 
	ord <- un_neighbor(un_matrix,500)
	ord_cor <- cor(ord,use="pairwise.complete.obs")
	image(ord_cor)

}
