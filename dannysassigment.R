prime<-function(number){
	for(x in 2:sqrt(number)){
		if(number%%x==0)stop("Not a prime number. It divides by: ",x)
	}
	cat("It's prime as hell! \n")
}

rowmeans<-function(mat){
	for(r in 1:nrow(mat)){
			avg<-0
			for(n in 1:ncol(mat)){
				avg<-avg+(mat[r,n]/ncol(mat))
			}
			cat("Mean of values in row ",r," is ",avg,"\n")
	}
}

hypotenuse<-function(a,b){
	c<-sqrt(a*a+b*b)
	print(c)
}

led<-function(on_seq){
	screen<-matrix("-",10,10)
	for(i in on_seq){
		r<-i%/%10+1
		k<-i%%10
		screen[r,k]<-i
	}
	print(screen)
}

faculty<-function(number){
	res<-1
	for(i in 1:number){
		res<-res*i
	}
	cat(res,"\n")
}

save_image_from_url <- function(image_url,image_file){
	con <- file(image_file,"wb")
	writeBin(getBinaryURL(image_url),con,useBytes=TRUE)
	close(con)
}


r_is_made_for_searching_google <- function(searchterm){
	#setwd(paste("D:/My_searches/",searchterm,sep=""))
	website <- getURL(paste("http://www.google.nl/images?q=",searchterm,sep=""))
	cat("Website is ",nchar(website),"characters at begining.\n")
	results <- regexpr("http(.+?).jpg",website)
	cat("First match at ",results[[1]],"\n")
	i<-1
	while(results>0){
		if(attr(results,"match.length")<256){
			cat(i,"OK!",substr(website,results[[1]],results[[1]]+attr(results,"match.length")-1),"\n")
			save_image_from_url(substr(website,results[[1]],results[[1]]+attr(results,"match.length")-1),paste("image",i,".jpg",sep=""))
			website <- substr(website,results[[1]]+attr(results,"match.length"),nchar(website))
		}else{
			website <- substr(website,results[[1]]+4,nchar(website))
		}
		i<-i+1
		results <- regexpr("http(.+?).jpg",website)
	}
}