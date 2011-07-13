#copy-paste into R for analysis
workflow.appriopriateMarkers <- function(){
	setwd("D:/data")
	library(basicQtl)
	library(qtl)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- toGenotypes(expressionMatrix,margin=0.5,genotypes=c(1,0),overlapInd=0, verb=T)
	recoMatrix <- recombinationCount(brassica_genotypes,flip=0)
	recoMatrixFlipped <- recombinationCount(brassica_genotypes,flip=1)
	res <- rep(0,length(rownames(recoMatrix)))
	names(res) <-  rownames(recoMatrix)
	for(x in 1:nrow(recoMatrix)){
		i <- recoMatrix[x,]
		j <- recoMatrixFlipped[x,]
		res[names(which(j[which(i>24)]<24))] <- res[names(which(j[which(i>24)]<24))]+1
	}
	hist(res)
	genotypes[which(res>780&&res<880),] <- 1-genotypes[which(res>780&&res<880),]
}

workflow.parental <- function(){
	setwd("D:/data/parental")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2)
	crossSimulated <- toGenotypes(ril,use="simulated",verbose=TRUE,debugMode=1,treshold=0.05,margin=15)
	
	design <- sort(ril$rils$phenotypes[1,],decreasing=TRUE)
	weight <- c(rep(0,73),rep(0,73))
	modellikelihood(design[c(1:10,35000:35010),],rep(3,146),c(rep(0,73),rep(0,73)))
	
	setwd("D:/GenesForSeedQualityWageningen")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2,sep="\t")
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2,groupLabels=c(0,0,0,1,1,1))
	ril2 <- ril
	ril2$rils$phenotypes <- ril2$rils$phenotypes[,-1]
	ril2$rils$phenotypes <- matrix(as.numeric(ril2$rils$phenotypes),nrow(ril$rils$phenotypes),(ncol(ril$rils$phenotypes)-1))
	ril2$parental$phenotypes <- ril2$parental$phenotypes[,-1]
	ril2$parental$phenotypes <- matrix(as.numeric(ril2$parental$phenotypes),nrow(ril$parental$phenotypes),(ncol(ril$parental$phenotypes)-1))
	rownames(ril2$rils$phenotypes)<- 1:nrow(ril2$rils$phenotypes)
	colnames(ril2$rils$phenotypes)<- 1:ncol(ril2$rils$phenotypes)
	rownames(ril2$parental$phenotypes)<- 1:nrow(ril2$parental$phenotypes)
	colnames(ril2$parental$phenotypes)<- 1:ncol(ril2$parental$phenotypes)
	cross <- toGenotypes(ril2,use="simulated",treshold=0.05,overlapInd = 0, proportion = c(50,50), margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	
	setwd("D:/data/from slave/data")
	ril <- readFiles(verbose=TRUE,debugMode=2,sep="")
	ril <- preprocessData(ril,verbose=TRUE,debugMode=2,groupLabels=c(0,0,1,1))
	cross <- toGenotypes(ril,use="simulated",treshold=0.05,overlapInd = 0, proportion = c(50,50), margin=15, minChrLength=5, verbose=TRUE, debugMode=1, max.rf=0.26, min.lod=0)
	
exp_p <- read.table(file="children_phenotypes_rp.txt" ,sep="\t")
exp_p <- read.table(file="parental_phenotypes_rp.txt" ,sep="\t")
	
	
	setwd("D:/data/bc")
	library(pheno2geno)
	ril <- readFiles(verbose=TRUE,debugMode=2)
	ril <- preprocessData(ril,groupLabels = c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1),verbose=TRUE,debugMode=2)
	crossMean <- toGenotypes(ril,use="simulated",treshold=0.01,overlapInd = 0, proportion = c(25,50,25), minChrLength=5, verbose=TRUE, debugMode=1)
	
	png("tolerateordereddefault.png", width=1200, height=1200)
	plot.rf(formLinkageGroups(cross,reorgMarkers=TRUE))
	dev.off()
 }

ril <- convertToGenotypes.internal(ril, 0, 0.01, c(50,50), 15, TRUE, 1)






setwd("D:/data/wilco data")
parental <- read.table(file="mapping_probes.txt",sep="\t",header=T,row.names=1)
children <- read.table(file="children.txt",sep="\t",header=T,row.names=1)
ril <- NULL
ril$rils$phenotypes <- children[,20:120]
ril$parental$phenotypes <- parental[,5:14]
ril <- preprocessData(ril, c(1,1,1,1,1,0,0,0,0,0))
cross <- toGenotypes(ril, use="simulated", splitMethod="EM", treshold=0.001, verbose=TRUE,debug=1)
chil2 <- parchil

setwd("D:/data/wilco data")
children <- read.table(file="children.txt",sep="\t",header=T,row.names=1)
parchil <- children[,8:120]
colnames(parchil) <- cnames[,2]
parchil <- parchil[,-which(colnames(parchil)=="Pimp_6_1")]
parchil <- parchil[,-which(colnames(parchil)=="MM_6_2")]
parchil <- parchil[,-which(colnames(parchil)=="RIL_278_6")]
tom <- tom[,-which(colnames(tom)=="..bot1703.CEL")]
tom <- tom[,-which(colnames(tom)=="..bot1707.CEL")]
tom <- tom[,-which(colnames(tom)=="..bot2000.CEL")]
#parchil <- parchil[,-which(colnames(parchil)=="RIL_260_d")]

parchil[,which(colnames(parchil)=="RIL_292_6")], parchil[,which(colnames(parchil)=="RIL_308_d")]
colo <- matrix("grey",1,111)
colo[which(colnames(parchil)=="RIL_278_6")] <- "red"
colo[which(colnames(parchil)=="RIL_260_d")] <- "red"
colo[which(colnames(parchil)=="RIL_206_6")] <- "red"
group_6 <- grep("_6",colnames(parchil))
group_d <- grep("_d",colnames(parchil))
groups <- list(group_6,group_d)

setwd("D:/data/wilco data")
parchil <- read.table(file="expressions_log_norm_cor.txt",sep="\t",header=T,row.names=1)
population <- createPopulation(parchil[,11:108],parchil[,1:10])
population <- preprocessData(population, c(1,1,1,1,1,0,0,0,0,0))
cross <- toGenotypes(population, use="simulated", splitMethod="mean", treshold=0.01, verbose=TRUE,debug=1)

show.dist <- function(nr){
	print(colnames(genotypes)[nr])
	res <- apply(genotypes,2,cor,genotypes[,nr])
	hist(res)
	print(names(res)[which(res>0.4)])
}

	filename <- paste("tomato_exp_",i,".txt",sep="")
	cat("Processing:",filename,"\n")
	tom <- read.table(filename,sep="")
	tom <- tom[,c(-13,-17,-91)]
	j <- 1
	tom_ <- matrix(apply(tom,1,diffexp),nrow(tom),1)
	print(dim(tom_))
	rownames(tom_) <- rownames(tom)
	filename2 <- paste("tomato_exp_",i,"_corrected.txt",sep="")
	write.table(tom[names(sort(tom_[,1],decreasing=TRUE)[1:100000]),], file=filename2, sep="\t")
	cat("Done:",filename2,"\n")
	tom <- NULL
	tom_ <- NULL
	gc()
	gc()
	gc()
	gc()
	gc()

setwd("D:/GenesForSeedQualityWageningen/tomato")
tom_c <- read.table("batchcorrected.txt",sep="\t")
tom <- read.table("200_000_markers_normalized.txt",sep="\t")
for(i in 1:12){
	filename <- paste("tomato_exp_",i,".txt",sep="")
	cat("Processing:",filename,"\n")
	tom <- read.table(filename,sep="",row.names=1)
	tom <- tom[,c(-13,-17,-91)]
	j <- 1
	tom_ <- matrix(apply(tom,1,diffexp),nrow(tom),1)
	print(dim(tom_))
	rownames(tom_) <- rownames(tom)
	filename2 <- paste("tomato_exp_",i,"_corrected.txt",sep="")
	write.table(tom[names(sort(tom_[,1],decreasing=TRUE)),], file=filename2, sep="\t")
	cat("Done:",filename2,"\n")
	tom <- NULL
	tom_ <- NULL
	gc()
	gc()
	gc()
	gc()
	gc()
}

for(i in 2:13){
	filename <- paste("tomato_exp_",i,"_corrected.txt",sep="")
	cat("Processing:",filename,"\n")
	tom <- read.table(filename,sep="\t")
	if(ncol(tom)==118){
		rownames(tom) <- tom[,1]
		tom <- tom[,-1]
	}
	colnames(tom) <- colnames(tom_c)
	cat(dim(tom_c),dim(tom),"\n")
	tom_c <- rbind(tom_c,tom)
	j <- 1
	tom_ <- matrix(apply(tom_c,1,diffexp),nrow(tom_c),1)
	print(dim(tom_))
	rownames(tom_) <- rownames(tom_c)
	tom_c <- tom_c[names(sort(tom_[,1],decreasing=TRUE)[1:250000]),]
	tom <- NULL
	tom_ <- NULL
	gc()
	gc()
	gc()
	gc()
	gc()
}
if(i%%5000==0)
norm <- NULL
for(i in 1:nrow(tom_c)){
	norm <- rbind(norm, normalizeQuantiles(tom_c[i,]))
	 cat(i,"/",nrow(tom_c),"\n")
}

for(i in 1:ncol(ara)){
	if(colnames(ara)[i]%in%anames[,1]){
		colnames(ara)[i] <- anames[which(anames[,1]==colnames(ara)[i]),2]
	}
}


cnames <- read.table("ril_labels_tomato.txt" ,sep="\t")
colnames(tom_c)<-as.character(cnames[which(cnames[,1]%in%colnames(tom_c)),2])
load("200_000_population.rd")
load()
tom_ <- matrix(unlist(tom_c),nrow(tom_c),ncol(tom_c))
rownames(tom_) <- rownames(tom_c)
colnames(tom_) <- colnames(tom_c)
tom_c <- tom_
tom_ <- NULL
gc()
gc()
gc()
gc()
population <- createPopulation(tom_c[,11:109],tom_c[,1:10])
population <- preprocessData(population, c(1,1,1,1,1,0,0,0,0,0))
cross <- toGenotypes(population, use="simulated", splitMethod="mean", treshold=0.01, verbose=TRUE,debug=2)
which(rownames(population$founders$phenotypes)[which(population$founders$RP$pval[1]<0.01)] %in% rownames(population$founders$phenotypes)[which(population$founders$RP$pval[2]<0.01)])

### arabidopsis

for(i in 1:19){
	filename <- paste("arabidopsis_exp_",i,".txt",sep="")
	cat("Processing:",filename,"\n")
	ara <- read.table(filename,sep="\t",header=F,row.names=1)
	colnames(ara) <- colnames(ara_c)
	cat(dim(ara_c),dim(ara),"\n")
	ara_c <- rbind(ara_c,ara)
	ara_ <- matrix(apply(ara_c,1,diffexp),nrow(ara_c),1)
	print(dim(ara_))
	rownames(ara_) <- rownames(ara_c)
	ara_c <- ara_c[names(sort(ara_[,1],decreasing=TRUE)[1:10000]),]
	ara <- NULL
	ara_ <- NULL
	gc()
	gc()
	gc()
	gc()
	gc()
}

diffexp <- function(cur_row){
	a <- mean(as.numeric(cur_row[17:24]),na.rm=TRUE)
	b <- mean(as.numeric(cur_row[189:196]),na.rm=TRUE)
	if(!(is.na(a))&&!(is.na(b))){
		invisible(abs(a-b))
	}else{
	invisible(0)
	}
}

for(i in 1:98){
	colnames(genotype)[i] <- rn[which(rn[,1]==i),3]
}

setwd("D:/GenesForSeedQualityWageningen/tomato")
require(pheno2geno)
load("200_000_population.rd")
map <- read.table("marker.txt",sep="\t",row.names=1)
map <- as.matrix(map)
genotype <- read.table("offspring_genotype.txt",sep="\t",header=T)
genotype <- as.matrix(genotype)
population <- intoPopulation(population, list(genotype,map), c("offspring$genotypes","maps$genetic"))
cross <- toGenotypes(population,genotype="real",orderUsing="maps_genetic")


setwd("D:/data/tomato")
require(pheno2geno)
population <- readFiles(verbose=T,debugMode=2)
population <- preprocessData(population,c(0,0,0,0,0,1,1,1,1,1))
cross <- toGenotypes(population,genotype="real",orderUsing="maps_genetic")
cross <- toGenotypes(ril,splitMethod="mean",genotype="simulated",minChrLength=0,treshold=0.5,margin=50,max.rf=10)
population_ <- removeIndividuals.internal(population,c("RIL_308_d","RIL_304_d","RIL_278_6"))
tom_c <- tom_c[,-which(colnames(tom_c)=="RIL_304_d")]

for(i in 1:ncol(tom_c)){
	colnames(tom_c)[i]<-cnames[which(cnames[,1]==colnames(tom_c)[i]),2]
}

result <- NULL
for(i in 1:nrow(pl)){
	print(sum(pl[,1]==pl[i,1]))
	if(sum(pl[,1]==pl[i,1])!=1){
		cur <- which(pl[,1]==pl[1,1])
		if(sum(pl[cur,3]==100)!=1){
			cur_ <- cur[which(pl[cur,3]==100)]
			if(sum(pl[cur_,4]==25)==1){
				result <- rbind(result,unlist(c(pl[cur_[which(pl[cur_,3]==100)],1],as.numeric(substr(pl[cur_[which(pl[cur_,3]==100)],2],9,10)),pl[cur_[which(pl[cur_,3]==100)],c(9,10)])))
			}
		}else{
			result <- rbind(result,unlist(c(pl[cur[which(pl[cur,3]==100)],1],as.numeric(substr(pl[cur[which(pl[cur,3]==100)],2],9,10)),pl[cur[which(pl[cur,3]==100)],c(9,10)])))
		}
		pl <- pl[-which(pl[,1]==pl[i,1]),]
	}else{
		result <- rbind(result,unlist(c(pl[i,1],as.numeric(substr(pl[i,2],9,10)),pl[i,c(9,10)])))
	}
	#result <- matrix(unlist(result),nrow(result),ncol(result))
}

sortTable <- function(table_,col_n){
	res_ <- NULL
	cur_ <- sort(table_[,col_n],decreasing=T)
	for(i in 1:nrow(table_)){
		res_ <- rbind(res_,table_[,which(table_[,col_n]==cur_[i])])
	}
	invisible(res_)
}