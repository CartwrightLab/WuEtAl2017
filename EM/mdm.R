# mixture of dirichlet-multinomials

# required libraries
library(mc2d)
library(MASS)

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

# fit a dirichlet-multinomial distribution
# using method of moments
mdmSingleMoments <- function(x,w=rep(1, nrow(x)),pseudoCounts=FALSE) {
	x <- as.matrix(x)
	if(pseudoCounts) {
		x <- x+1
	}
	# remove columns that have no size
	y <- colSums(x)
	oo <- (y > 0)
	x <- x[,oo]

	# remove rows that have no size
	n <- rowSums(x)
	x <- x[n>0,]
	n <- n[n>0]
	N <- nrow(x)
	# control for different sample sizes
	q <- x/n
	# calcuate proportions
	m <- cov.wt(q,w,method="ML")
	p <- m$center
	v <- m$cov
	# calcualte phi
	h <- weighted.mean(1/n,w)
	v <- v / (diag(p)-outer(p,p))
	phi <- (mean(v)-h)/(1.0-h)
	if(phi <= 0.0) {
		phi <- 1/(sum(w*n)+1)
	} else if( phi >= 1.0) {
		phi <- 1-1/(sum(w*n)+1)
	}
	pp <- vector("numeric",length(oo))
	pp[oo] <- p
	mdmParams(phi,pp)
}

# fit a mixture of dirichlet-multinomial distributions
# using method of moments
mdmMoments <- function(x,M=1L,w=rep(1, nrow(x)),pseudoCounts=FALSE,nstart=1) {
	if(!is.numeric(M) || M < 1) {
		stop("Parameter M must be a positive integer.")
	}
	M <- as.integer(M)

	# remove rows with no data
	x <- as.matrix(x)
	# apply the weights, assumes integer values
	x <- x[rep(1:nrow(x),w),]

	rownames(x) <- NULL
	n <- rowSums(x)
	x <- x[n>0,]
	n <- n[n>0]
	q <- x/n

	# calculate a kmeans clustering
	# this can vary between runs, so try 10 different starts
	ret <- kmeans(q,M,nstart=nstart)

	s <- split(data.frame(x),ret$cluster)
	f <- ret$size/nrow(x)
	params <- t(sapply(s,mdmSingleMoments))
	params[params < .Machine$double.eps/4] <- .Machine$double.eps/4
	names(f) <- paste("m",seq_len(M),sep="")
	rownames(params) <- names(f)
	colnames(params) <- c("phi",paste("p",seq_len(ncol(x)),sep=""))
	class(params) <- "mdmParams"
	list(f=f,params=params)
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

# remove empty columns, construct mask, and append column sizes
mdmAugmentData <- function(x,w=NULL) {
	x <- as.matrix(x)
	# remove empty columns
	y <- colSums(x)
	oo <- (y > 0)
	x <- x[,oo]
	n <- rowSums(x)
	r <- cbind(x,n,deparse.level=0)

	int <- do.call("interaction", c(unclass(as.data.frame(x)),drop=TRUE))
	y <- r[match(levels(int),int),]
	w <- mdmTabulateWeights(int,w)

	list(r=r,y=y,w=w,mask=oo)
}

# calculate the summary statistics
mdmSumStats <- function(x,w=NULL,augmented=FALSE) {
	if(augmented) {
		r <- list(r=x,mask=rep(TRUE,ncol(x)))
	} else {
		r <- mdmAugmentData(x)
	}
	mx <- max(r$r)
	u <- apply(r$r,2,function(y) mdmTabulateWeights(y,w,mx))
	s <- apply(u,2,function(y) rev(cumsum(rev(y))))
	return(list(s=s,mask=r$mask))
}

mdmSingleCore <- function(s,params,phTol,cycles,phAcc,traceLevel=0) {
	logTol <- log(10)*(phTol/-10.0)
	logAcc <- log(10)*(phAcc/-10.0)
	# if message is not specified create our own
	msg <- function(i,pars,ll) {
		ss <- sprintf("%0.6g", pars)
		ll <- sprintf("%0.16g", ll)
		z <- sprintf("    %d: %s (ll = %s)", i, paste(ss,collapse=" "),ll)
		message(z)
	}
	# use only the first row of params
	params <- as.vector(params)
	# setup variables
	N <- nrow(s)
	KK <- ncol(s)
	K <- KK-1
	nn <- sum(s[,KK])
	# stopping criteria of Narayanan (1991)
	#   NOTE: using observed information
	#   http://www.jstor.org/stable/2347605
	x2.stop <- qchisq(logTol,K,log.p=TRUE)
	x2.acc  <- qchisq(logAcc,K,log.p=TRUE,lower.tail=FALSE)

	if(traceLevel >= 2) message("Fitting dirichlet-multinomial:")

	# vectorization is easier if we transpose s
	ss <- s
	s <- t(s)
	j <- rep(0:(N-1), each=KK)

	# TODO: check to see if data fits phi=1

	ll2  <- -.Machine$double.xmax
	ll   <- -.Machine$double.xmax
	lln <- mdmSingleLogLikeCore(ss,params)

	nparams <- params
	for(cyc in 1:cycles) {
		params <- nparams
		# sometimes we need to correct for p being zero due to
		# numerical precision
		params[-1][params[-1] == 0.0] <- .Machine$double.eps/4

		ll2 <- ll
		ll <- lln
		phi <- params[1]
		p <- c(params[-1],1)

		# calculate gradient
		gf <- rowSums(s*(j-p)/(p+phi*(j-p)))
		gg <- rowSums(s/(p+phi*(j-p)))*(1.0-phi)
		gp <- gg[1:(K-1)]-gg[K]
		g <- c(sum(gf[-KK])-gf[KK], gp)

		# calculate hessian
		hh <- -rowSums(s/(p+phi*(j-p))^2)*(1.0-phi)^2
		h <- matrix(0,K,K)
		h[2:K,2:K] <-  (diag(hh[1:(K-1)],nrow=(K-1),ncol=(K-1))+hh[K])
		hh <- -rowSums(s*j/(p+phi*(j-p))^2)
		h[2:K,1] <- (hh[1:(K-1)]-hh[K])
		h[1,2:K] <- h[2:K,1]
		hh <- -rowSums(s*(j-p)^2/(p+phi*(j-p))^2)
		h[1,1] <- sum(hh[-KK])-hh[KK]

		# invert hessian, use pseudo-inverse incase it is singlular
		# TODO: determine if we can calculate the inverse directly
		# TODO: Woodbury matrix identity and block matrices
		ih <- ginv(-h)

		# calculate delta
		delta <- ih %*% g

		#delta <- solve(-h,g)

		# Test For Convergence
		# Sometimes when phi=0, the likelihood
		# is maximized, but g != 0.  Detect this and fix g.
		#
		gAdj <- g
		if(phi == 0.0 && g[1] < 0.0) {
			gAdj[1] <- 0.0
		}

		# this statistic is approximately chi-square distributed
		# so we can stop when it is statistically near enough to 0
		x2 <- t(gAdj) %*% ih %*% gAdj
		if(x2 >= 0 && x2 < x2.stop) {
			break
		}

		# calculate new parameters
		nparams <- params[1:K] + delta
		nparams <- c(nparams,1-sum(nparams[2:K]))
		lln <- mdmSingleLogLikeCore(ss,nparams)
		# if we overshoot, use an alternative search
		if(ll > lln) {
			# Use a simple fixed-point search to update p
			np <- p[-KK]*gg[-KK]
			np <- np/sum(np)

			# Approximate the likelihood using a generalized Newton
			# approximation:  f(phi) = -b Log[1 - phi] - (a + b) phi + k
			# where a, b, and k are set such that f(phi), f'(phi), and f''(phi)
			# all fit the log-likelihood surface.
			# As appropriate, f(phi) = -a Log[phi] - (a + b) (1 - phi) + k
			# is used instead.
			f1 <- g[1]
			f2 <- h[1,1]
			fpp <- f2*(1-phi)*phi

			if((f2 >= 0.0 && fpp >= f1) || (f2 <= 0.0 && fpp <= f1)) {
				b <- f2*(1-phi)^2
				a <- -f1+fpp
			} else {
				a <- f2*phi^2
				b <- f1+fpp
			}
			if(a+b != 0.0) {
				# will maximize the approximate function as long as f2 is negative
				nphi <- a/(a+b)
			} else {
				nphi <- phi
			}

			nparams <- c(nphi, np)

			lln <- mdmSingleLogLikeCore(ss,nparams)
			if(ll > lln) {
				# For badly fitting data, the above approximation can fail.
				# Use a slower and safer fixed-point method instead.
				# Try a simple Aitken accelleration procedure
				z <- rowSums(s)
				u <- z*phi-phi^2*gf

				nphi0 <- phi*u[KK]/(sum(u[-KK])-phi*(sum(u[-KK])-u[KK]))
				gf <- rowSums(s*(j-p)/(p+nphi0*(j-p)))
				u <- z*nphi0-nphi0^2*gf
				nphi1 <- nphi0*u[KK]/(sum(u[-KK])-nphi0*(sum(u[-KK])-u[KK]))
				gf <- rowSums(s*(j-p)/(p+nphi1*(j-p)))
				u <- z*nphi1-nphi1^2*gf
				nphi2 <- nphi1*u[KK]/(sum(u[-KK])-nphi1*(sum(u[-KK])-u[KK]))

				nphi <- nphi0-(nphi1-nphi0)^2/((nphi2-nphi1)-(nphi1-nphi0))

				nparams <- c(nphi, np)
				lln <- mdmSingleLogLikeCore(ss,nparams)
				if(ll > lln) {
					# If Aitken accelleration has failed resort to using
					# the first iteration of the fixed-point method.
					nparams <- c(nphi0, np)
					lln <- mdmSingleLogLikeCore(ss,nparams)
				}
			}
		}
		#llnull <- sum(s[-KK,]*log(p))
		if(traceLevel >= 2) msg(cyc-1,params,ll)
	}
	if(traceLevel >= 2) msg(cyc,params,ll)
	list(ll=ll,params=params,obsInformation=-h,covar=ih,score=g,cycles=cyc)
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

mdmLogLikeCore <- function(r,w,f,params) {
	# setup variables
	N <- nrow(r)
	KK <- ncol(r)
	n <- r[,KK]
	K <- KK-1
	M <- length(f)

	# sanity checks
	if(KK != ncol(params)) {
		stop("ncol(params) != ncol(r).");
	}
	if(M != nrow(params)) {
		stop("nrow(params) != length(f).");
	}

	# calculate component probabilities for each row
	pr <- matrix(0,N,M)
	mx <- max(n)
	j <- rep(seq.int(0,mx-1),each=KK)
	for(m in 1:M) {
		phi <- params[m,1]
		p <- c(params[m,-1],1)
		if(phi == 1.0) {
			g <- max.col(r,ties.method="first")
			p <- c(log(f[m])+log(p),-1/0)
			pr[,m] <- p[g]
		} else {
			cache <- cbind(0,matrix(log(p+phi*(j-p)),nrow=KK))
			cache <- apply(cache,1,cumsum)
			ww <- sapply(1:KK,function(x) cache[1+r[,x],x])
			pr[,m] <- log(f[m])+rowSums(ww[,-KK])-ww[,KK]
		}
	}

	ww <- log(rowSums(exp(pr)))
	pr <- pr-ww
	list(ll = sum(w*ww), w = exp(pr))
}


# fit a dirichlet-multinomial distribution
# using maximum likelihood
#TODO: paramter sanity check
mdmSingle <- function(x,phi=NULL,p=NULL,phTol=100,phAcc=40,cycles=1000,traceLevel=0) {
	s <- mdmSumStats(x)
	mask <- c(TRUE,s$mask)
	if(is.null(phi)) {
		# use the method of moments to start the search
		iniParams <- mdmSingleMoments(x)
	} else {
		iniParams <- mdmParams(phi,p)
	}
	params <- iniParams[1,mask]
	if(!all(0.0 <= params & params <= 1.0 )) {
		warning("invalid initial parameter string")
		return(NULL)
	}

	ret <- mdmSingleCore(s$s,params,phTol,cycles,phAcc,traceLevel)
	KK <- length(ret$params)
	K <- KK-1

	L <- length(mask)
	retParams <- vector("numeric",L)
	names(retParams) <- colnames(iniParams)
	retParams[mask] <- ret$params

	mask2 <- mask
	mask2[max(which(mask))] <- FALSE

	freeParams <- retParams[mask2]

	covar <- matrix(0,L,L,dimnames=list(names(retParams),names(retParams)))
	oo <- which(mask)
	pos <- cbind(rep(oo,L),rep(oo,each=L))
	A <- diag(rep(1,KK))
	A[KK,2:K] <- -1
	cv <- cbind(rbind(ret$covar,0),0)
	oo <- which(mask)
	L <- length(oo)
	pos <- cbind(rep(oo,L),rep(oo,each=L))
	covar[pos] <- (A %*% cv %*% t(A))

	info <- ret$obsInformation
	score <- ret$score
	names(score) <- names(freeParams)
	colnames(info) <- names(freeParams)
	rownames(info) <- names(freeParams)

	vecParams <- retParams
	#retParams <- t(retParams)
	#rownames(retParams) <- "m1"

	list(ll=ret$ll,
		params=vecParams,freeParams=freeParams,
		covar=covar,obsInformation=info,
		score=score,cycles=ret$cycles)
}

# fit a mixture of dirichlet-multinomial distributions
# using maximum likelihood
mdm <- function(x,f=1L,phi=NULL,p=NULL,phTol=100,cycles=1000,
	cyclesInner=NULL,phAcc=40,traceLevel=0L,fixStart=FALSE) {

	if(missing(x)) {
		stop("x not specified.")
	}

	if(length(f) == 1 && f == 1) {
		ret <- mdmSingle(x,phi,p,phTol=phTol,
				cycles=cycles,traceLevel=traceLevel)
		return(ret)
	}
	# calculate tolerance
	logTol <- (phTol/-10.0)*log(10)
	logAcc <- (phAcc/-10.0)*log(10)
	# augment data rows
	r <- mdmAugmentData(x)
	mask <- c(TRUE,r$mask)
	mx <- max(r$y)
	if(is.null(phi)) {
		# use the method of moments to start the search
		M <- if(length(f) == 1) f else length(f)
		mm <- mdmMoments(x,M)
		f <- mm$f
		iniParams <- mm$params
		params <- iniParams
	} else {
		iniParams <- mdmParams(phi,p)
		params <- iniParams[,mask]
	}
	N <- sum(r$w)
	NN <- nrow(r$y)
	KK <- ncol(r$y)
	K <- KK-1
	M <- length(f)
	ndf <- M*KK-1

	if(M < 2) {
		stop("f must contain two or more values.")
	}
	f <- f/sum(f)
	iniF <- f

	x2.stop <- qchisq(logTol, ndf,log.p=TRUE)
	x2.acc  <- qchisq(logAcc, ndf,log.p=TRUE,lower.tail=FALSE)

	#calculate cycles for single-model estimate
	if(is.null(cyclesInner)) {
		cyclesInner <- max(KK,cycles/M)
	}

	#setup likelihood
	ll   <- -.Machine$double.xmax
	ll2  <- -.Machine$double.xmax
	ll3  <- -.Machine$double.xmax
	llStar  <- -.Machine$double.xmax
	llStar2 <- -.Machine$double.xmax

	# a bad starting point can cause all sorts of problems
	# use initial parameters to define components
	# follow with method of moments
	if(fixStart) {
            w <- mdmLogLikeCore(r$y,r$w,f,params)
            wsum <- colSums(r$w*w$w)
            for(m in 1:M) {
                params[m,] <- mdmSingleMoments(r$y[,-KK],r$w*w$w[,m])
            }
            f <- wsum/N
	}

	# calculate likelihood and weights
	w <- mdmLogLikeCore(r$y,r$w,f,params)
	newParams <- params
	newParamsA <- params
	newF <- f
	newFA <- f
	for(cyc in seq_len(cycles)) {
		params <- newParams
		f <- newF

		# sometimes we need to correct for p and f being zero due to
		# numerical precision
		params[,-1][params[,-1] == 0.0] <- .Machine$double.eps/4
		f[f == 0.0] <- .Machine$double.eps/4

		ll3 <- ll2
		ll2 <- ll
		ll <- w$ll

		if(traceLevel>0)
			message(sprintf("\n**** Cycle %d ****",cyc))

		if(traceLevel>0) {
			if(identical(cyc,1L)) {
				message(sprintf("ll = %0.16g", ll))
			} else {
				message(sprintf("ll = %0.16g (%+0.16g)", ll, ll-ll2))
			}
			cat("Parameters:\n")
			print(cbind(f,params))
		}
		if(ll2 > ll && (traceLevel>0 || !isTRUE(all.equal(ll2,ll)))) {
			warning(sprintf("Cycle %d: Log-Likelihood decreased by %0.16g!",
				cyc, ll2-ll))
		}

		# E-step
		# calculate weights
		wsum <- colSums(r$w*w$w)

		#### Calculate the Score Vector for each row
		So <- matrix(0,NN,ndf)
		# f params
		for(m in seq_len(M-1)) {
			So[,m] <- w$w[,m]/f[m]-w$w[,M]/f[M]
		}
		# p and phi params
		jKK <- rep(0:(mx-1),each=KK)
		jK <- rep(0:(mx-1),each=K)
		for(m in seq_len(M)) {
			pos <- M+K*(m-1)
			phi <- params[m,1]

			# phi
			p <- c(params[m,-1],1)
			aj <- (jKK-p)/(p+phi*(jKK-p))
			aj <- matrix(aj,nrow=KK)
			aj <- cbind(0,aj)
			scoreCache <- apply(aj, 1, cumsum)
			gg <- matrix(0,NN,KK)
			for(k in seq_len(KK)) {
				gg[,k] <- scoreCache[1+r$y[,k],k]
			}
			So[,pos] <- w$w[,m]*(rowSums(gg[,-KK])-gg[,KK])

			# p
			p <- c(params[m,-1])
			aj <- (1.0-phi)/(p+phi*(jK-p))
			aj <- matrix(aj,nrow=K)
			aj <- cbind(0,aj)
			scoreCache <- apply(aj, 1, cumsum)
			gg <- matrix(0,NN,K)
			for(k in seq_len(K)) {
				gg[,k] <- scoreCache[1+r$y[,k],k]
			}
			So[,pos+(1:(K-1))] <- w$w[,m]*(gg[,-K]-gg[,K])
		}
		wSo <- r$w*So

		#### Caclulate the Observed Information Matrix
		#### Calculate the Observed Full Data Information Matrix
		Io <- matrix(0,ndf,ndf)
		Bo <- Io
		# f-component of the matrix
		Bo[1:(M-1),1:(M-1)] <- wsum[M]/(f[M]*f[M])
		for(m in seq_len(M-1)) {
			Bo[m,m] <- Bo[m,m] + wsum[m]/(f[m]*f[m])
		}
		# phi and p compoents
		for(m in seq_len(M)) {
			pos <- M+K*(m-1)
			phi <- params[m,1]
			p <- params[m,-1]
			s <- mdmSumStats(r$y,r$w*w$w[,m],TRUE)
			#TODO: optimize out the creaion of ss?
			ss <- t(s$s[,-KK])
			# p x p
			hh <- -rowSums(ss/(p+phi*(jK-p))^2)*(1.0-phi)^2
			Bo[pos+(1:(K-1)),pos+(1:(K-1))] <- -hh[K]
			for(k in 1:(K-1)) {
				Bo[pos+k,pos+k] <- -hh[k]-hh[K]
			}
			# phi x p
			hh <- -rowSums(ss*jK/(p+phi*(jK-p))^2)
			Bo[pos,pos+(1:(K-1))] <- -hh[-K]+hh[K]
			Bo[pos+(1:(K-1)),pos] <- Bo[pos,pos+(1:(K-1))]
			# phi x phi
			#TODO: optimize out the creaion of ss?
			ss <- t(s$s)
			p <- c(params[m,-1],1)
			hh <- -rowSums(ss*(jKK-p)^2/(p+phi*(jKK-p))^2)
			Bo[pos,pos] <- -sum(hh[-KK])+hh[KK]

			for(ki in 0:(K-1)) {
				for(kj in 0:(K-1)) {
					Io[pos+ki,pos+kj] <- Bo[pos+ki,pos+kj] -
						sum(wSo[,pos+ki]*So[,pos+kj]/w$w[,m])
				}
			}
		}
		# f x (phi,p) components
		for(m in seq_len(M-1)) {
			posI <- m
			for(k in seq_len(K)) {
				posJ <- M+K*(m-1)+(k-1)
				Io[posI,posJ] <- -sum(wSo[,posJ])/f[m]
			}
			for(k in seq_len(K)) {
				posJ <- M+K*(M-1)+(k-1)
				Io[posI,posJ] <- sum(wSo[,posJ])/f[M]
			}
		}

		# E(i)*E(j) and make symmetrical matrix
		for(posI in seq_len(M*KK-1)) {
			Io[posI,posI] <- Io[posI,posI] + sum(wSo[,posI]*So[,posI])
			posJ <- posI+1
			while(posJ < M*KK) {
				Io[posI,posJ] <- Io[posJ,posI] <-
					Io[posI,posJ] + sum(wSo[,posI]*So[,posJ])
				posJ <- posJ+1
			}
		}
		# Invert the observed information matrix and test for convergence.
		# Use pseudo-inverse incase Io is (near)-singular.  It seems to
		# work well.
		#Vo <- try(ginv(Io))
		Vo <- ginv(Io)
		if(inherits(Vo, "try-error")) {
			x2 <- -1
		} else {
			# Sometimes phi=0 maximizes the likelihod without gradient being
			# zero.  Detect and fix.
			Se <- colSums(wSo)
			gAdj <- Se
			mZero <- which(params[,1] == 0.0 & Se[M+K*seq(0,M-1)] < 0.0)
			gAdj[M+K*(mZero-1)] <- 0.0

			#
			x2 <- t(gAdj) %*% Vo %*% gAdj
		}
		if(traceLevel>0) {
			message(sprintf("Convergence:  x2 = %0.16g", x2))
		}

		# stopping criteria of Narayanan (1991)
		#   http://www.jstor.org/stable/2347605
		if((x2 >= 0 && x2 < x2.stop) || ll == ll2) {
			break
		}

		# Aitken stopping criteria of McLachlan and Krishnan (2008) 4.9
		if(FALSE && cyc > 2) {
		  llStar2 <- llStar
		  ci <- (ll-ll2)/(ll2-ll3)
		  llStar <- ll2 + (ll-ll2)/(1.0-ci)
		  if(traceLevel>0)
		  	cat(sprintf("ll* = %0.16g\n",llStar))
		  if(abs(llStar2 - llStar)/N < exp(logTol)) {
		  	break
		  }
		}

		# Optimize parameters for each component
		for(m in seq_len(M)) {
			# E-Step
			s <- mdmSumStats(r$y,r$w*w$w[,m],TRUE)
			# M-Step
			ret <- mdmSingleCore(s$s,params[m,],phTol,cyclesInner,phAcc,traceLevel)
			newParams[m,] <- ret$params
		}
		# M-Step
		newF <- wsum/N

		# Accelerate?
		if(x2 >= 0 && x2 < x2.acc) {
			np <- c(newF[-M],t(newParams[,-KK]))
			op <- c(f[-M],t(params[,-KK]))
			fp <- op + (Vo %*% Bo) %*% (np-op)

			newFA <- fp[seq_len(M-1)]
			newFA <- c(newFA,1-sum(newFA))
			newParamsA[,-KK] <- matrix(fp[-(1:(M-1))],nrow=M,ncol=K,byrow=TRUE)
			newParamsA[,KK] <- 1.0-rowSums(newParamsA[,2:K,drop=FALSE])

			# Control Overshooting
			if(any(newFA <= 0.0 | 1.0 <= newFA )) {
				newFA <- newF
			}

			for(m in seq_len(M)) {
				if(any(newParamsA[m,] < 0.0 | 1.0 < newParamsA[m,])) {
					newParamsA[m,] <- newParams[m,]
				}
			}
			w <- mdmLogLikeCore(r$y,r$w,newFA,newParamsA)
			if(w$ll > ll) {
				newF <- newFA
				newParams <- newParamsA
			} else {
				w <- mdmLogLikeCore(r$y,r$w,newF,newParams)
			}
		} else {
			w <- mdmLogLikeCore(r$y,r$w,newF,newParams)
		}
	}
	# calculate a mask for free parameters
	L <- length(mask)
	mask2 <- mask
	mask2[max(which(mask))] <- FALSE
	freeMask <- c(rep(TRUE,M-1),FALSE,rep(mask2,M))
	varMask <- c(rep(TRUE,M),rep(mask,M))
	# apply names to estimated parameters
	retF <- f
	names(retF) <- paste("m",seq_len(M),sep="")
	retParams <- matrix(0,M,L,dimnames=list(names(retF),colnames(iniParams)))
	retParams[,mask] <- params

	# construct vectors of parameters and free parameters
	vecParams <- c(f,t(retParams))
	names(vecParams) <- c(
		paste("f", seq_len(M),sep="."),
		paste(colnames(retParams),rep(1:M,each=ncol(retParams)),sep=".")
	)
	freeParams <- vecParams[freeMask]

	# construct transition matrix for coverting freeParam covar matrix
	U <- matrix(0,nrow=(M+M*(K+1)),ncol=(1+M*K+M-1))
	for(m in seq_len(M-1)) {
		U[m,m] <- 1
		U[M,m] <- -1
	}
	for(m in seq_len(M)) {
		for(k in seq_len(K)) {
			U[M+(K+1)*(m-1)+k,M-1+(K)*(m-1)+k] <- 1
		}
		for(k in 2:K) {
			U[M+(K+1)*(m-1)+K+1,M-1+(K)*(m-1)+k] <- -1
		}
	}
	U[c(M,M+(K+1)*(1:M)),1+M*K+M-1] <- 1

	# calculate covariance matrix for vecParams
	cv <- cbind(rbind(Vo,0),0)
	cv <- (U %*% cv %*% t(U))
	covar <- matrix(0,length(vecParams),length(vecParams),dimnames=list(names(vecParams),names(vecParams)))
	oo <- which(varMask)
	L <- length(oo)
	pos <- cbind(rep(oo,L),rep(oo,each=L))
	covar[pos] <- cv

	colnames(Io) <- names(freeParams)
	rownames(Io) <- names(freeParams)

	list(ll=ll,f=retF,params=retParams,
		vecParams=vecParams,freeParams=freeParams,
		obsInformation=Io,covar=covar, score=Se,
		cycles=cyc
	)
}

# Paul et al (2005) 10.1002/bimj.200410103
mdmSingleInfo <- function(x,phi,p=NULL) {
	params <- mdmParams(phi,p)
	# calculate alphas
	alphas <- mdmAlphas(params,total=TRUE)

	x <- as.matrix(x)
	N <- nrow(x)
	K <- ncol(x)
	KK <- K+1

	A <- alphas[KK]
	phi <- params[1]

	n <- rowSums(x)
	nt <- tabulate(n)
	m <- length(nt)
	w <- which(nt>0)

	s <- matrix(0,K,m)
	for(k in seq.int(K)) {
		a <- alphas[k]
		b <- alphas[KK]-a
		alog <- cumsum(log(a+seq.int(0,m-1)))
		blog <- rev(c(0,cumsum(log(b+seq.int(0,m-2)))))
		for(i in w) {
			lch <- lchoose(i,seq.int(1,i))
			glog <- lch + alog[seq.int(1,i)] + blog[seq.int(m-i+1,m)] -
				(lgamma(a+b+i)-lgamma(a+b))
			ss <- nt[i]*rev(cumsum(rev(exp(glog))))
			s[k,seq_along(ss)] <- s[k,seq_along(ss)] + ss
		}
	}

	j <- rep(seq.int(0,m-1),each=K)
	p <- alphas[-KK]/A

	hm <- rowSums(s/(j+alphas[-KK])^2)
	h <- matrix(0,K,K)
	hh <- A*A*hm
	h[2:K,2:K] <- (diag(hh[1:(K-1)],nrow=(K-1),ncol=(K-1))+hh[K])

	hh <- rowSums(s*j/(j+alphas[-KK])^2)*(A+1)^2
	h[2:K,1] <- hh[1:(K-1)]-hh[K]
	h[1,2:K] <- h[2:K,1]

	hh <- p*rowSums(s*(2*j+alphas[-KK]-p)/(alphas[-KK]+j)^2)
	j <- seq.int(0,ncol(s)-1)
	s <- rev(cumsum(rev(nt)))
	hz <- sum(s*(2*j+alphas[KK]-1)/(alphas[KK]+j)^2)
	h[1,1] <- -(sum(hh)-hz)*(A+1)^3
	h
}
