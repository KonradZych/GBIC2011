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