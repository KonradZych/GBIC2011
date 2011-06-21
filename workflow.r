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
cnames <- read.table("ril_labels_tomato.txt" ,sep="\t")
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

diffexp <- function(cur_row){
	a <- mean(as.numeric(cur_row[8:12]),na.rm=TRUE)
	b <- mean(as.numeric(cur_row[13:17]),na.rm=TRUE)
	if(!(is.na(a))&&!(is.na(b))){
		invisible(abs(a-b))
	}else{
	invisible(0)
	}
}

load("batchcor.rd")
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
cross <- toGenotypes(population, use="simulated", splitMethod="mean", treshold=0.2, verbose=TRUE,debug=1)
which(rownames(population$founders$phenotypes)[which(population$founders$RP$pval[1]<0.01)] %in% rownames(population$founders$phenotypes)[which(population$founders$RP$pval[2]<0.01)])
