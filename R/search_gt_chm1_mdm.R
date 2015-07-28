#! /usr/bin/Rscript --vanilla
# filetype=R

# Rscript search-gt-mdm.R ${numTrial} ${numCom}${isDirty} $file 
traceLevel <- 0
upperLimit <- 150
lowerLimit <- 10
#set.seed(17761980)
getScriptPath <- function(){

    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.dir <- dirname(script.name)

    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}


script.dir = getScriptPath()
suppressPackageStartupMessages(source( file.path(script.dir,"mdm.R") )) 
options(warn=1)

fitmdm.mle.many <- function(x,n,k) {
	fv <- rdirichlet(n,c(rep(5/(k-1),k-1),95))
	am <- matrix(c(1,1,0.1),1,3)
	for(i in 2:k) {
		j <- 1+0.5*trunc((i-1)/2)
		if(i %% 2 == 0) {
			am <- rbind(am,c(j,1,0.1))
		} else {
			am <- rbind(am,c(1,j,0.1))
		}
	}
	fiv <- rdirichlet(k*n,10*am)
	
	res <- list()
	m <- -.Machine$double.xmax
	for(i in 1:n) {
		f <- fv[i,]
		od <- sort(0.2*runif(k)^2)
		a <- fiv[seq(k*(i-1)+1,k*i),,drop=FALSE]*(1-od)/od

		a <- a[order(rowSums(a)),]
		
		if(k == 1) {
			f <- 1
		}
		
		u <- try(mdm(x,f,a,phTol=150,cycles=1000,cyclesInner=200,
				phAcc=40,traceLevel=traceLevel,fixStart=FALSE))
		if(inherits(u, "try-error")) {
			message(sprintf("  Trial %d: FAILED (max ll = %0.16g)", i,m))
		} else {
			m <- max(m,u$ll)
			message(sprintf("  Trial %d: ll = %0.16g (max ll = %0.16g)", i, u$ll,m))
			res <- c(res,list(u))
		}
	}
	res
}

#########################################################


pid <- Sys.getpid()
args <- commandArgs(trailingOnly=TRUE)
dirty_data <- FALSE

if(length(args) < 3) {
	message("  Usage: search-gt-mdm (num of trials) (num of compoents[d]) (data file).")
	message("  To do a search using dirty data, append 'd' to the number of components.")
	quit(status=1)
}

nn <- as.numeric(args[1])
k <- args[2]
dirty_data <- grepl("^\\d+[Dd]$",k)
if(dirty_data) {
	k <- as.numeric(sub("[Dd]$","",k))
} else {
	k <- as.numeric(k)
}

dat <- read.delim(args[3],header=TRUE)

sub_name <- gsub("base_counts_(.*_CHM1.*)_byref_.*txt","\\1" ,args[3])
if(sub_name == args[3]){
	message("Can't find the pattern: ", sub_name)
	sub_name <- ""
}


x <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
x <- data.matrix(x)
row.names(x) <- dat$pos
x <- x[dat$callby == 1 & ((dat$snp == 1 & dat$snpdif == 0) | dirty_data), ]
n <- rowSums(x)
oo <- lowerLimit <= n & n <= upperLimit
x <- x[oo,]
n <- n[oo]

message(sprintf("Searching for maximum likelihood of %d-component model....",k))
message(sprintf("Using%s data with %dx to %dx coverage....",
	(if(dirty_data) " dirty" else ""), lowerLimit, upperLimit))

ofile <- sprintf(  
	"gt_mdm_%s_%d%s_results_%d.RData",sub_name, k, if(dirty_data) "D" else "", pid)

res <- fitmdm.mle.many(x,nn,k)
ll <- sapply(res,function(x) x$ll)

resName <- sprintf("res_%d", pid)
resNameK <- sprintf((if(dirty_data) "resK_%dd" else 
	"resK_%d"),k,pid)

res <- res[rev(order(ll))]

assign(resName,res)
assign(resNameK,resName)

cat("Saving results to", ofile, "...\n")

save(list=c("resName",resName,resNameK),file=ofile)

