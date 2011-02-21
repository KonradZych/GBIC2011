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

#clean - Removing NA from data matrix by replacing it with 0
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

#makebinary - making a binary matrix out of matrix containing data e.g. from gene expression data -> gene expressed/not expressed
#matrix_to_be_made_binary: matrix containing data
#returns: output - matrix containing 0s (value below the treshold) and 1 (value above the treshold) currently, the treshold is median
makebinary<-function(matrix_to_be_made_binary){
	output<-matrix(0,nrow(matrix_to_be_made_binary),ncol(matrix_to_be_made_binary))
	#resulting matrix should have the same labels as input
	rownames(output)<-rownames(matrix_to_be_made_binary, do.NULL = FALSE)
	colnames(output)<-colnames(matrix_to_be_made_binary, do.NULL = FALSE)
	#using median as a treshold value
	for(i in 1:ncol(matrix_to_be_made_binary)){
	    tres=median(matrix_to_be_made_binary[,i])
		output[,i]<-(matrix_to_be_made_binary[,i]>tres)
	}
	output
}

#qtlbyttest - Basic single marker mapping by using a t.test statistic (For RIL)
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
		if(mean(pheno_class1)>=mean(pheno_class2)){
		output <- c(output,-log10(t.test(pheno_class1,pheno_class2)[[3]]))	
		}else{
		output <- c(output,log10(t.test(pheno_class1,pheno_class2)[[3]]))
		}		
	}
	output
}


#heatmapqtl - Creates data file concerning all traits and also image of it
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

#pathway - creates probable pathway using values for specified marker
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

#makecorvector - double mode:
#line 100!
#0 - using one matrix, producing matrix of correlations between neighbor values in one direction 
#1 - using two matrices, one as above, second to add exact values
#input - one or two matrices of data, two_matrices (number 0 or 1)
#function returns matrix suitable for make_topo_pallete
makecorvector <- function(first_matrix,two_matrices=0,second_matrix=NULL){
	first_matrix_cor <- cor(first_matrix,use="pairwise.complete.obs")
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

#make_topo_pallete - making nice collor pallete using RGB function, values >0 are red and less transparent the higher their are, <0 blue and less transparent the lower their are
#input - matrix produced by makecorvector
#returns - matrix of colors in RGB,alpha format, max value for color and alpha =255
make_topo_pallete <- function(result_matrix){
	cur_mean <- mean(abs(result_matrix))
	cur_sd <- mean(sd(abs(result_matrix)))
	cur_range <- abs(max(abs(result_matrix))-min(abs(result_matrix)))
	topo_pallete<-matrix(0,nrow(result_matrix),ncol(result_matrix))
	for(i in 1:nrow(result_matrix)){
		for(j in 1:ncol(result_matrix)){
			if(result_matrix[i,j]>=0){
				if(result_matrix[i,j]<(cur_mean-cur_sd)){
					topo_pallete[i,j]<-rgb(0,0,abs(result_matrix[i,j]/cur_range)*255,255-abs(result_matrix[i,j]/cur_range)*255,maxColorValue=255)
				}else if((result_matrix[i,j]>(cur_mean+cur_sd))){
					topo_pallete[i,j]<-rgb(abs(result_matrix[i,j]/cur_range)*255,0,0,abs(result_matrix[i,j]/cur_range)*255,maxColorValue=255)
				}else{
					topo_pallete[i,j]<-rgb(0,abs(result_matrix[i,j]/cur_range)*255,0,55-abs(result_matrix[i,j]/cur_range)*55,maxColorValue=255)
				}	
			}else{
				cur<-abs(result_matrix[i,j])
				if(cur>(cur_mean+cur_sd)){
					topo_pallete[i,j]<-rgb(0,0,abs(cur/cur_range)*255,abs(cur/cur_range)*255,maxColorValue=255)
				}else if((cur<(cur_mean-cur_sd))){
					topo_pallete[i,j]<-rgb(abs(cur/cur_range)*255,0,0,255-abs(cur/cur_range)*255,maxColorValue=255)
				}else{
					topo_pallete[i,j]<-rgb(0,abs(cur/cur_range)*255,0,55-abs(cur/cur_range)*55,maxColorValue=255)
				}	
			}
		
		}
	}
	topo_pallete
}

#un_intonumeric - function formating data to be usable, specificly for one file, but can be easily adapted to other, mind the comment below in the text
#input - matrix of data
#returns matrix with specified characters
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

#un_recombination - instead of cor use recombination, very very slow (triple for - Danny won't like it)
#input - matrix of data
#return - matrix of recombination values between COLUMNS
un_recombination<-function(chrom_matrix){
	output <- matrix(0,ncol(chrom_matrix),ncol(chrom_matrix))
	#triple for, jupi!
	for(i in 1:ncol(chrom_matrix)){
		cat("Analysing column ",i,"\n")
		for(j in 1:ncol(chrom_matrix)){
			for(k in 1:nrow(chrom_matrix)){
				if(is.na(chrom_matrix[i,j])||is.na(chrom_matrix[i,j])){
					rec <- 0
				}else if(chrom_matrix[k,j]==chrom_matrix[k,i]){
					rec <- 0
				}else{
					rec <- 1
				}
				output[i,j] <- output[i,j] + rec
			}
		}	
	}
	#beacuase we want to be able to use the same functions as for corelaction, recombination values must be scaled
	output <- (100-output)/100
#line 200!
	output
}

#un_row_score - counting score for specified vector, which means, how many values w are above specified treshold, divided by it's length
#input - vector of values, treshold (number 0-1)
#return - row score (number)
un_row_score <- function(cor_matrix_row,treshold=0.7){
	row_score <- 0
	for(j in 1:length(cor_matrix_row)){
		if(cor_matrix_row[j]>treshold){
			row_score <- row_score + 1/length(cor_matrix_row)
		}
	}
	row_score
}

#un_drop_markers - removing columns, that are higly correlated with more that specified percentage of others
#input - matrix of data, treshold (number 0-1)
#return - matrix of data of the same type
un_drop_markers <- function(chrom_matrix,treshold=0.25){
	cor_matrix <- un_recombination(chrom_matrix)
	result <- apply(cor_matrix,1,un_row_score)
	i<-1
	while(max(result)>treshold){
		chrom_matrix <- chrom_matrix[,-(which(result==max(result)))]
		cor_matrix <- un_recombination(chrom_matrix)
		result <- apply(cor_matrix,1,un_row_score)
		print(i)
		print(max(result))
		i<-i+1
	}
	chrom_matrix
}

#un_best_clustering - needs improvment really bad! quadro-for!:D - making spceified number of clustering of data
#then producing a matrix of points, udes later for further classification
#input - matrix of data, number of iterations(int), number of groups(int)
#return - matrix of numbers (0-nr_iterations)
un_best_clustering <- function(chrom_matrix,nr_iterations,groups=10){
	cor_matrix <- un_recombination(chrom_matrix)
	print("un_best_clustering starting")
	res <- NULL
	map <- matrix(0,nrow(cor_matrix),ncol(cor_matrix))
	print("iteration starting")
	#clustering with k-means
	for(i in 1:nr_iterations){
		r <- kmeans(cor_matrix,groups)
		res <- rbind(res,(as.numeric(r[[1]])))
	}
	print("iteration done, starting pointing system")
	#matrix of points
	for(i in 1:nr_iterations){
		for(j in 1:groups){
			for(k in which(res[i,]==j)){
				for(l in which(res[i,]==j)){
					map[k,l] <- map[k,l] + 1
				}
			}
		}
	}
	print("pointing done, returning output")
	#matrix should inherit colnames from input
	colnames(map)<-colnames(chrom_matrix, do.NULL = FALSE)
	map
}

#un_order_chromosome - ordering markers inside one group (chromosome)
#input - matrix of data (specified fragment to be sorted inside)
#return - names of columns in sorted order
un_order_chromosome <- function(chrom_matrix){
	cat(ncol(chrom_matrix)," markers\n")
	output<-chrom_matrix
	#sorting is made in number of iterations equal to number of columns
	for(i in 1:ncol(chrom_matrix)){
		cat("Starting iteration ",i,"\n")
		chrom_cor_matrix <- un_recombination(output)
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
		output <- chrom_matrix[,result]
	}
	print("Iterations done,saving result")
	output <- colnames(chrom_matrix[,result])
	output	
}

#un_neighbor - heart of analysis!
#input: matrix of data with wrongly ordered columns, to be clustered, sorted inside groups, nr_iterations (int) groups(int)
#line 300!
#return: matrix of the same data with rigth order of columns
un_neighbor <- function(un_matrix,nr_iterations=1000,groups=5){
	print("Running un_best_clustering")
	r <- un_best_clustering(un_matrix,nr_iterations,groups)
		for(i in 1:nrow(r)){
		for(j in 1:ncol(r)){
			if(r[i,j]<20){r[i,j]<-0}
		}
	}
	print("un_best_clustering done, running un_best_clustering")
	r <- un_best_clustering(r,nr_iterations,groups)
	print("un_best_clustering done, final clustering")
	r<-kmeans(r,groups)
	print("Final clustering done, running un_order_chromosome")
	res <- NULL
	for(i in 1:groups){
		cat("Segregating chromosome: ",i,"\n")
		cur <- un_order_chromosome(un_matrix[,which(r[[1]]==i)])
		res <- cbind(res,un_matrix[,cur])
	}
	print("un_order_chromosome done, returning result")
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

#examples to be used in package
#
#1 - basic qtl map, using data from gene expression (phenotypes.txt) and genotyping (genotypes.txt)  
#
makebinary_test()
setwd("G:/git")
phenotypes <- as.matrix(read.table("phenotypes.txt", sep=""))
genotypes <- as.matrix(read.table("genotypes.txt", sep=""))
result_binary<-heatmapqtl(makebinary(clean(phenotypes)),genotypes)
covariance_matrix <- makecorvector(genotypes,1,result_binary)
color_pallete <- make_topo_pallete(covariance_matrix)
persp(result_binary,col=color_pallete)
#obviously, one can skip this fancy coloring, but, hell, it's nice looking
#
#2 - obtaining simple pathawy from observing marker 100 peak
#
setwd("G:/git")
phenotypes <- as.matrix(read.table("phenotypes.txt", sep=""))
genotypes <- as.matrix(read.table("genotypes.txt", sep=""))
path<-pathway(result_binary, 100)
#
#3 - recreating genemap from messed data
#
setwd("G:/git")
un_matrix <- un_intonumeric(as.matrix(read.table("genotypes_multitrait.txt", sep="", header=TRUE)))
un_result<-un_drop_markers(un_matrix)
ord <- un_neighbor(un_result,1000,5)
ord_recombination <- un_recombination(ord)
image(ord_cor)