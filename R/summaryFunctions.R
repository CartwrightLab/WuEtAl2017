require(mc2d)

extractMaxModel<- function(path, isCEU=True){
        if(isCEU){
            listRData<- list.files(path=path, pattern="hets.+RData") 
        } else {
           listRData<- list.files(path=path, pattern="homos_alt.+RData") 
        }
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



calculateEachLikelihood<- function(maxModel, fullData, lowerLimit, upperLimit, numData=NULL, isCEU=TRUE){

    whichIsDirty <- grepl("_[0-9]D",names(maxModel))
    dataRef<- parseData(fullData, lowerLimit, upperLimit, dirtyData=FALSE, isCEU)
    dataRefDirty<- parseData(fullData, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

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


parseData<- function(dat, lowerLimit, upperLimit, dirtyData, isCEU=TRUE){
    dataRef <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
    dataRef <- data.matrix(dataRef)
    row.names(dataRef) <- dat$pos
    if(isCEU){
        dataRef <- dataRef[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]
    } else {
        dataRef <- dataRef[dat$callby == 1 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]
    }
    
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
	

plotqq<- function(z, ff, outerText){
    mains = c("Reference Allele", "Alternate Allele", "Error")

    for(i in 1:3) {
        qqplot(z[,i],ff[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
        abline(0,1)
    #   print(ad.stat.k(ff[,i],z[,i]))
    }
    mtext(outerText, side=3, outer=T, line=-2, cex=2)

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


## copied from mdm.R

# generate a random sample of DM observations
# n = number of observations
# m = a vector (or scalar) of observation sizes
# phi = a vector (or scalar) of dispersion parameters
# p = a matrix (or vector) of proportions
#   if p is NULL, phi is assumed to contain alpha (scale) parameters
rdm <- function(n, m, phi, p=NULL) {
	params <- mdmParams(phi,p)
	if(any(params[,1] < 0.0 | 1.0 < params[,1])) {
		stop("phi must be in [0,1].")
	}
	if(any(params[,-1] <= 0.0 | 1.0 <= params[,-1])) {
		stop("p must be in (0,1).")
	}
	params[params[,1] < .Machine$double.eps/2,1] <- .Machine$double.eps/2
	p <- params[,-1,drop=FALSE]
	
	alphas <- mdmAlphas(params)
	# choose initial 
	y <- rmultinomial(n,1,p)
	# update params, conditional on what has occurred
	ny <- nrow(y)
	na <- nrow(alphas)
	if(na != ny) {
		n1 <- ny %/% na
		n2 <- ny %% na
		u <- rep(seq_len(na),n1)
		if(n2 > 0) {
			u <- c(u,seq_len(n2))
		}
		a <- y + alphas[u,]
	} else {
		a <- y + alphas
	}
	# choose following
	y+rmultinomial(n,m-1,rdirichlet(n,a))	
}

# generate a random sample of mixture of DM distributions
# n = number of observations
# m = a vector (or scalar) of observation sizes
# f = the mixture proportions
# phi = a vector (or scalar) of dispersion parameters
# p = a matrix (or vector) of proportions
#   if p is NULL, phi is assumed to contain alpha (scale) parameters
rmdm <- function(n, m, f, phi, p=NULL) {
	params <- mdmParams(phi,p)
	k <- nrow(params)
	if(length(f) != k) {
		stop("The length of 'f' and number of rows in params must be equal.")
	}
	
	# generate the mixture
	q <- rmultinomial(1, n, f)
	mix <- rep.int(seq_len(k), q)
	mix <- sample(mix)
	# generate the samples
	x <- rdm(n,m,params[mix,])
	# use the rownames to store the mixture categories
	rownames(x) <- mix
	x
}

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

# convert parameters to a parameter vector
# phi = a vector (or scalar) of dispersion parameters
# p = a matrix (or vector) of proportions
#   if p is NULL, phi is assumed to contain alpha (scale) parameters
mdmParams <- function(phi, p=NULL) {
	if(inherits(phi, "mdmParams")) {
		return(phi)
	}
	if(!is.null(p)) {
		if(is.vector(p)) {
			p <- t(p)
		}
		p <- p/rowSums(p)
	} else {
		if(is.vector(phi)) {
			a <- t(phi)
		} else {
			a <- phi
		}
		A <- rowSums(a)
		phi <- 1.0/(A+1.0)
		p <- a/A		
	}

	colnames(p) <- paste("p", seq_len(ncol(p)),sep="")
	v <- cbind(phi,p)
	class(v) <- "mdmParams"
	v
}

`[.mdmParams` <- function(x, i, j, ...) {
  y <- NextMethod(.Generic)
  class(y) <- .Class
  y
}

