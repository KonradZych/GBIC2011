#WORK here
check_parameters <- function(functionName, numericParameters, booleanParameters, debugModeParameters, verbose=FALSE, debugMode=0){
	if(verbose || debugMode==1 || debugMode==2 || debugMode==3) cat("Checking parameters for:",functionName,"\n")
	if(numericParameters){if(!lapply(numericParameters,is.numeric))stop("One of parameters you provided to",functionName,"is not numeric. Check help file.\n")}
	if(booleanParameters){if(!lapply(booleanParameters,))stop("Verbose you provide to",functionName,"must be boolean. Check help file.\n")}
	if(debugModeParameters){if(!lapply(debugModeParameters,))stop("DebugMode you provide to",functionName,"must be 0,1,2 or 3. Check help file.\n")}
}

isBoolean <- function(x){
	if(x!=TRUE&&x!=FALSE){ return(FALSE)}
	else{ return(TRUE) }
}

isDebug <- function(x){
	if(!pmatch(x,c(0,1,2,3))){ return(FALSE)}
	else{ return(TRUE) }
}

flipMarkers <- function(genotypicMatrix){
	brassicaReco <- recombinationCount(genotypicMatrix)
	brassicaRecoflipped <- recombinationCount(genotypicMatrix,flip=1)
	brassicaRecorows <- apply(brassicaReco,1,mean)
	brassicaRecoflippedrows <- apply(brassicaRecoflipped,1,mean)
	result <- NULL
	for(i in 1:length(brassicaRecoflippedrows)){
		cat(i,"\n")
		if(brassicaRecorows[i]<brassicaRecoflippedrows[i]){
			result <- rbind(result, genotypicMatrix[i,])
		}else{
			result <- rbind(result, (1-genotypicMatrix)[i,])
		}
	}
	result
}