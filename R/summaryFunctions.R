

extractMaxModel<- function(path){

	listRData<- list.files(path=path, pattern="hets.+RData") 
	maxModel<- list()

	for(i in 1:length(listRData)){

		load(paste(path, listRData[i], sep=""))
		data<- get(resName)
		if(length(data) != 0){
			index <- which.max(sapply(data, function(x) {x$ll} ))
			maxModel[[i]]<- data[[index]]
			if (is.null(maxModel[[i]][["f"]]) ){ # Fix 1 parameter model
				maxModel[[i]][["f"]] <- 1
				maxModel[[i]]$params <- t(as.matrix(maxModel[[i]]$params))
			}
		}
		else{
			maxModel[[i]]<- "NULL"
		}
		
		subName <- gsub("gt_mdm_(hets_.*)_result.*\\.RData","\\1" ,listRData[i])
		names(maxModel)[[i]]<- subName
	}
	return(maxModel)
}



calculateEachLikelihood<- function(maxModel, fullData, lowerLimit, upperLimit, numData=NULL){

    whichIsDirty <- grepl("_[0-9]D",names(maxModel))
    dataRef<- parseData(fullData, lowerLimit, upperLimit, dirtyData)
    dataRefDirty<- parseData(fullData, lowerLimit, upperLimit, dirtyData=TRUE)

    maxLikelihood <- vector(length=length(maxModel), mode="list")
    names(maxLikelihood)<- names(maxModel)
    
    for( d in 1:length(maxModel)){

        if(whichIsDirty[d]){
            data<- dataRefDirty
        }
        else{
            data<- dataRef
        }
        maxLikelihood[[d]]<- calculateEachLikelihoodOneModel(maxModel[[d]], data)
        
#         params<- maxModel[[d]]$param
#         if(!is.matrix(params)){
#             params<- matrix(params, nrow=1)
#         }
#         numParams <- nrow(params)
#         numData<- nrow(data)
        
#         maxLikelihood[[d]]<- matrix(nrow=numData, ncol=numParams)
#         for (i in 1:numData){
#             x <- matrix(data[i,], nrow=1)
#             sum_stat <- mdmSumStatsSingle(x)
#             stat <- sum_stat$s
# 
#             for(p in 1:numParams){
#                 maxLikelihood[[d]][i,p]<- mdmSingleLogLikeCore(stat, params[p,])
#             }
#         }

    }
    
    return(maxLikelihood)
}

calculateEachLikelihoodOneModel<- function(model, data){
    numData<- nrow(data)
    
    params<- model$param
#         if(!is.matrix(params)){
#             params<- matrix(params, nrow=1)
#         }
    numParams <- nrow(params)
    
    maxLikelihood<- matrix(nrow=numData, ncol=numParams)
    for (i in 1:numData){
        x <- matrix(data[i,], nrow=1)
        sum_stat <- mdmSumStatsSingle(x)
        stat <- sum_stat$s

        for(p in 1:numParams){
            maxLikelihood[i,p]<- mdmSingleLogLikeCore(stat, params[p,])
        }
    }
    return(maxLikelihood)
}


parseData<- function(dat, lowerLimit, upperLimit, dirtyData){
    dataRef <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
    dataRef <- data.matrix(dataRef)
    row.names(dataRef) <- dat$pos
    dataRef <- dataRef[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]
    n <- rowSums(dataRef)
    oo <- lowerLimit <= n & n <= upperLimit
    dataRef <- dataRef[oo,]
    n <- n[oo]
    return(dataRef)
}


getMaxLikelihoodProp<- function(maxLikelihoodList){
	result<- sapply(maxLikelihoodList, function(x){
                prop <- table(apply(x,1,which.max))
                return(prop/sum(prop))
			})
	return(result)
}
    

parseDataFP<- function(dat, lowerLimit, upperLimit){
	dataRef <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
	dataRef <- data.matrix(dataRef)
	row.names(dataRef) <- dat$pos
	dataRef <- dataRef[dat$callby == 2 & dat$snp == 0, ]
	n <- rowSums(dataRef)
	oo <- lowerLimit <= n & n <= upperLimit
	dataRef <- dataRef[oo,]
	n <- n[oo]
	return(dataRef)
}
	
parseDataIndex<- function(dat, index, lowerLimit, upperLimit){
    dataRef <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
    dataRef <- data.matrix(dataRef)
    row.names(dataRef) <- dat$pos
    dataRef <- dataRef[index, ]
    n <- rowSums(dataRef)
    oo <- lowerLimit <= n & n <= upperLimit
    dataRef <- dataRef[oo,]
    n <- n[oo]
    return(dataRef)
}
	
	
######################################################################
##### Modified from EM script
######################################################################
	

# convert a parameter vector to alphas
# v = a paramter vector contain phi and p
mdmAlphas <- function(v,total=FALSE) {
	if(is.vector(v)) {
		v <- t(v)
	}
	phi <- v[,1]
	p <- v[,-1,drop=FALSE]
	at <- ((1.0-phi)/phi)
	a <- p * at
	colnames(a) <- paste("a", seq_len(ncol(a)),sep="")
	if(total) {
		a <- cbind(a,aa=at)
	}
	rownames(a) <- NULL
	a
}


mdmAugmentDataSingle <- function(x,w=NULL) {
# 	x <- as.matrix(x)
	# remove empty columns
#  	y <- colSums(x)
#  	oo <- (y > 0)
#  	x <- as.matrix(x[,oo])
	n <- rowSums(x)
	r <- cbind(x,n,deparse.level=0)
	
	int <- do.call("interaction", c(unclass(as.data.frame(x)),drop=TRUE))
	y <- r[match(levels(int),int),]
	w <- mdmTabulateWeights(int,w)
	
	list(r=r,y=y,w=w)#,mask=oo)
}

# calculate the summary statistics
mdmSumStatsSingle <- function(x,w=NULL,augmented=FALSE) {
	r <- mdmAugmentDataSingle(x)
	mx <- max(r$r)
	u <- apply(r$r,2,function(y) mdmTabulateWeights(y,w,mx))
	s <- apply(u,2,function(y) rev(cumsum(rev(y))))
	return(list(s=s,mask=r$mask))
}


# tabulate weights
mdmTabulateWeights <- function(bin,w=NULL,nbins= max(1L, bin, na.rm = TRUE)) {
    if (!is.numeric(bin) && !is.factor(bin))
        stop("'bin' must be numeric or a factor")
    if (typeof(bin) != "integer") 
        bin <- as.integer(bin)
    if (nbins > .Machine$integer.max) 
        stop("attempt to make a table with >= 2^31 elements")
    nbins <- as.integer(nbins)
    if (is.na(nbins)) 
        stop("invalid value of 'nbins'")
	if(is.null(w)) {
		u <- .Internal(tabulate(bin, nbins))
	} else {
		u <- sapply(split(w,factor(unclass(bin),levels=1:nbins)),sum)
		names(u) <- NULL
	}
	u
}


mdmSingleLogLikeCore <- function(s,params) {
	# setup variables
	N <- nrow(s)
	KK <- ncol(s)
	n <- s[1,KK]
	K <- KK-1
	if(any(is.nan(params)) || any(params < 0.0 | 1.0 < params)) {
		return(-1.0/0.0)
	}
	if(KK != length(params)) {
		stop("ncol(s) != length(params)")
	}
	# vectorization is easier if we transpose s
	s <- t(s)
	s[KK,] <- -s[KK,]
	p <- c(params[-1],1)
	phi <- params[1]
	if(phi == 1.0) {
		tol <- 16*.Machine$double.eps
		if(!isTRUE(all.equal(0,sum(s[,1],tolerance=tol)))) {
			# if this is true, then the likelihood is -infinity
			return(-1/0)
		}
		return(sum(s[,1]*log(p)))
	}
	j <- rep(seq.int(0,N-1), each=KK)
	ll <- sum(s*log(p+phi*(j-p)))
	if(is.nan(ll)) {
		return(-1/0)
	}
	return(ll)
}
