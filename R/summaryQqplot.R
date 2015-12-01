#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

isCEU <- FALSE
isCEU <- TRUE


dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10

latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

if(isCEU){
    pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
    subNameList<- c(
    "CEU10_C10", "CEU10_C21",
    "CEU11_C10", "CEU11_C21",
    "CEU12_C10", "CEU12_C21",
    "CEU13_C10", "CEU13_C21"
    )

    fullTitleList<- c(
    "CEU 2010 Chromosome 10", "CEU 2010 Chromosome 21",
    "CEU 2011 Chromosome 10", "CEU 2011 Chromosome 21",
    "CEU 2012 Chromosome 10", "CEU 2012 Chromosome 21",
    "CEU 2013 Chromosome 10", "CEU 2013 Chromosome 21"
    )
} else {
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
}
setwd(pwd)


p<- 8
p<- 1
# for (p in 2:4){
# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# setwd(fullPath)
hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

rowSumDataRef<- rowSums(dataRef)
colSumDataRef<- colSums(dataRef)
rowSumDataRefDirty<- rowSums(dataRefDirty)

freqDataRef<- dataRef/rowSumDataRef
freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty

subName2<- gsub("_C", " Chr", subName)

maxModel<- extractMaxModel(fullPath)
whichIsDirty <- grepl("_[0-9]D",names(maxModel))
header<- gsub("hets_" , "", names(maxModel) )
header<- gsub(".*_", "", header)
header<- gsub("$", " components", header)
header<-paste(subName2, header)

qqplotFile<- file.path(latexDir, paste0("qqPlots_", subName, ".pdf") )


pdf(file=qqplotFile, width=12, height=6, title=qqplotFile)
par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), 
    cex.main=1.2^4,cex.lab=1.2^2, 
    omi=c(0,0,0.5,0) )

# mains = c("Reference Allele", "Alternate Allele", "Error")

#simple multinomial (no bias)
p1<- (colSumDataRef[1]+colSumDataRef[2])/2
prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
ll <- sum(log(prob)*colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, paste0("\nMultinomial ", subName2, "\n") )


#multinomial (with ref bias)
prob <- colSumDataRef / sum(colSumDataRef)
ll <- sum(log(prob)* colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, paste0("\nBiased Multinomial ", subName2, " \n"))

for( m in 1:length(maxModel)) {

    #dirichlet-multinomial mixture
    mModel<- maxModel[[m]]
    if(! whichIsDirty[m]){
        b <- rmdm(length(rowSumDataRef), rowSumDataRef, mModel$f,
                    phi=mModel$params[,1], p=mModel$params[,2:4])
        z <- b/rowSums(b)
        ff<- freqDataRef    
    } else{
#         b <- rmdm(length(rowSumDataRefDirty), rowSumDataRefDirty, mModel$f,
#                     phi=mModel$params[,1], p=mModel$params[,2:4])
#         ff<- freqDataRefDirty
        next
    }
    main<- paste0("\nMixture of Dirichlet Multinomial ", header[m] ," \n")
    if(m==1){
        main<- gsub(" 1 components", "", main)
        main<- gsub("Mixture of ", "", main)
    }
    plotqq(z, ff, main)

}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")
# } # match (p in 1:length(subNameList) ){


########################################
# originally run it with dm.R
# b <- rmdm(length(n),n,m3$param.p,m3$param.a)
# #NOTE: difference between $param $parmas
# # mdmParams should be $params but $params don't have class
# class(m3$params)<- "mdmParams" ## Wont work if class not forwords
# b <- rmdm(length(n),n,m3$f,m3$params)

empiricalQQ<- function(data){

    rowSum<- rowSums(data)
    colSum<- colSums(data)


    weight<- table(rowSum) / sum(table(rowSum))

    for(c in 10:150){
    
# #simple multinomial (no bias)
# p1<- (colSumDataRef[1]+colSumDataRef[2])/2
# prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
# ll <- sum(log(prob)*colSumDataRef)
# # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
# 
# b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
# z <- b/rowSums(b)
# 

dmultinomial(matrix(c(0:10,10:0,rep(0,11)), ncol=3), 10, prob)

dmultinomial( matrix(c(rep(0,11),10:0,0:10  ), ncol=3), 10, prob)


#collapse to binomial




    }

}

dmdm<- function(data, params){

    nCat<- NROW(params)
    alpha<- params[,1]
 
    apply(data,1,function(x){
        x
    
    })

    
    
}

allCombNG2<- function(N, mModel){
    site<- c(0,0)
    for(i in 0:N){
        M<- N-i
#         for (j in 0:M){
#             for(k in 0:(M-j)){
                j<- N - i
                cat(i," ",j," ", "\n")
                site<- rbind(site, c(i,j) )
#                 cat(dMDM(site,mModel) , "\n")
#             }
#         }
    }
    return(site[-1,])

}


allCombN<- function(N){
    site<- c(0,0,0)
    for(i in 0:N){
        M<- N-i
        for (j in 0:M){
                k<- N - i -j
#                 cat(i," ",j," ",k, "\n")
                site<- rbind(site, c(i,j,k) )

        }
    }
    return(site[-1,])

}
allX<-allCombN(2)
allX<-allCombN(5)
allX2<-allCombNG2(5, mModel)

mModel<- maxModel[[m]]

mModel$f <- c(1)
mModel$params<- rbind( c(0.001,0.5,0.5 ))


mModel$f <- c(0.5,0.5)
mModel$params<- rbind( c(0.001,0.5,0.4,0.1) , c(0.005,0.3,0.3,0.4))

params<-mModel$params
f<- mModel$f



dMDM<- function(x, mModel){
    
    alpha<- mdmAlphas(mModel$params)
    f<- mModel$f
    N<- sum(x)

    nCat<- NROW(alpha)
    
    rAll<- vector(length=nCat)
    for(cIndex in 1:nCat){
        sumAlpha<- sum(alpha[cIndex,])
        norm<- lgamma(sumAlpha) - lgamma(sumAlpha+N)
        r<- norm
        for(i in 1:3){
            if(x[i] > 0){
                for( read in 0:(x[i]-1)){
#                      cat(i, "\t", read, "\n")
                    r <- r + log(alpha[cIndex,i]+read)
                }
            }
        }
        
#         for( read in 0:(N-1)){
# #          cat(i, "\tN: ", read, "\n")
#                 r <- r - log(sumAlpha+read)
#         }
#         cat(prop[cIndex], "\t", r, "\n")
        rAll[cIndex] <- log(f[cIndex]) + r
    }
    return(rAll)
}


dMDM2<- function(x, mModel){
    
    alpha<- mdmAlphas(mModel$params)
    f<- mModel$f
    N<- sum(x)
    nCat<- NROW(alpha)
    
    rAll<- vector(length=nCat)
    for(cIndex in 1:nCat){
        sumAlpha<- sum(alpha[cIndex,])
        norm<- lgamma(sumAlpha) - lgamma(sumAlpha+N)
        r<- norm
        
        for(i in 1:3){
            r<- r + lgamma(alpha[cIndex,i]+x[i]) - lgamma(alpha[cIndex,i])
        }
#         cat(prop[cIndex], "\t", r, "\n")
        rAll[cIndex] <- log(f[cIndex]) + r
    }
    rAll <- rAll + lfactorial(N) - sum(lfactorial(x))
    return(sum(exp(rAll)))
}


dMDM3<- function(x, mModel){
    
    if(!is.matrix(mModel$params[,-1])){
#         mModel$params<- t(as.matrix(mModel$params))
        params<- t(as.matrix(mModel$params[,-1]))
    } else {
        params<- mModel$params[,-1]
    }
    allPhi<- mModel$params[,1]
    f<- mModel$f
    N<- sum(x)
    nCat<- length(allPhi)
    nGroup<- NCOL(params)
    
    rAll<- vector(length=nCat)
    for(cIndex in 1:nCat){
        phi<- allPhi[cIndex]
        oneMPhi<- 1-phi
#         log( 1 / prod( (1-phi+phi*(0:(N-1))) )   )
        norm <- -sum( log (oneMPhi+phi*(0:(N-1)))  ) 
#         * prod( 0.2*(1-phi)+ phi*(0:(k-1)) ) * prod ( 0.8*(1- phi) + phi*(0:(N-k-1))) 
        r<- norm
        for(i in 1:nGroup){
            if(x[i] > 0){
#                 cat(cIndex, " " , i , " " , params, "\n")
                r<- r + sum(log( params[cIndex,i]*oneMPhi + phi*(0:(x[i]-1)) ) )
            }
        }
#         cat(prop[cIndex], "\t", r, "\n")
        rAll[cIndex] <- log(f[cIndex]) + r
    }
    rAll <- rAll + lfactorial(N) - sum(lfactorial(x))
    return(sum(exp(rAll)))
}

mModel<- maxModel[[m]]

dmultinomial(allX,prob=c(0.5,0.4,0.1)) + dmultinomial(allX,prob=c(0.3,0.3,0.4))
dmultinomial(allX,prob=c(0.4996293251, 0.4996293251, 0.0007413498) )


N<- 150
allX<-allCombN(N)
dMDM(allX[2,], mModel)
dMDM2(allX[2,], mModel)
dMDM3(allX[2,], mModel)


prob<- apply(allX,1,function(y){dMDM3(y,mModel=mModel)})
prob
sum(prob)


collapseProb3<- function(allX, prob){
    cbind((collapseProb(allX, prob, 1))[,1],
        (collapseProb(allX, prob, 2))[,1], 
        (collapseProb(allX, prob, 3)) )
}

collapseProb<- function(allX, prob, group){
    max <- max(allX)
    mProb<- vector(length= max+1)
    for(i in 0:max){
        mProb[i+1]<- sum(prob[allX[,group]==i])
    }
    result<- cbind(cumsum(mProb), seq(0,1,length=max+1) )
    return(result)
}

collapseProbV2<- function(allX, prob, group){
    max <- max(allX)
#     mProb<- vector(length= max+1)
#     for(i in 0:max){
#         mProb[i+1]<- sum(prob[allX[,group]==i])
#     }
    
    mProb<- sapply(0:max, function(x){
        sum(prob[allX[,group]==x])
    })
    
    result<- cbind(cumsum(mProb), seq(0,1,length=max+1) )
    return(result)
}



# collapseProb3(allX, prob)

# d<-max(allX)+1
# probArray<- array(-1,c(d,d,d))
# dd<- allX+1
# for(i in 1:NROW(allX)){
#     print(allX[i,])
#         probArray[dd[i,1],dd[i,2],dd[i,3]] <- normProb[i]
# }


mModel<- maxModel[[m]]
mModel$f <- c(0.3,0.7)
mModel$params<- rbind( c(0.1,0.3,0.6,0.1) , c(0.05,0.8,0.1,0.1) )

mModel$f <- c(0.5,0.5)
mModel$params<- rbind( c(0.00, 0.4996293251, 0.4996293251, 0.0007413498),
        c(0.00, 0.4996293251, 0.4996293251, 0.0007413498) )

mModel$f <- c(1)
mModel$params<- c(0.00, 0.4996293251, 0.4996293251, 0.0007413498)
        
        
mModel$f <- c(0.5,0.5)
mModel$params<- rbind( c(0.001,0.5,0.4,0.1) , c(0.005,0.3,0.3,0.4))



N<- 10
allX<- allCombN(N)
prob<- apply(allX,1,function(y){dMDM3(y,mModel=mModel)})
prob2<- apply(allX,1,function(y){dmultinomial(y,prob=c(0.49,0.49,0.02))})
collapseProb3(allX, prob)


readWeightTemp<- table(rowSumDataRef) / sum(table(rowSumDataRef))
readWeight<- cbind(as.numeric(names(readWeightTemp)), readWeightTemp)
# readWeight<- readWeight[1:50, ]
# readWeight[,2] <- readWeight[,2] / sum(readWeight[,2] )

allCumProb<- vector(length=150, mode="list")

#simple multinomial (no bias)
p1<- (colSumDataRef[1]+colSumDataRef[2])/2
multiProb <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)

for( N in readWeight[,1]){
#     cat(N, "\n")
    allX<- allCombN(N)
    prob<- apply(allX,1,function(y){dMDM3(y,mModel=mModel)})
#    prob<- apply(allX,1,function(y){dmultinomial(y,prob=multiProb )})
    allCumProb[[N]]<- collapseProb3(allX, prob)
}

fetchIndex<- 1
size<- NROW(freqDataRef)
size<- 5000
plotData<- cbind(freqDataRef[1:size, fetchIndex],rep(-1,size))
site<- 1

for(f in 1:size){
    N<- rowSumDataRef[f]
    freq<- freqDataRef[f,fetchIndex]
    expect<- 0

    for( i in 1:NROW(readWeight)){
        cumProb<- allCumProb[[ readWeight[i,1] ]]
        index<- findInterval(freq, cumProb[,4])
        cp<- cumProb[index, fetchIndex] 
        expect<- expect + readWeight[i,2]*cp
    }


    # cumProb<- allCumProb[[ N ]]
    # index<- findInterval(freq, cumProb[,4])
    # cp<- cumProb[index, fetchIndex]
    # expect <- cp


    plotData[f, 2]<- expect


}

plot(plotData, xlim=c(0,1), ylim=c(0,1))

ac<- cbind(sort(plotData[,1]), cumsum(sort(plotData[,1]))/sum(plotData[,1]) )
# aaa<- vector(length=NROW(plotData))
# for(i in 1:NROW(plotData)){
#     index<- findInterval(plotData[i,1], ac[,1])
#     aaa[i]<- ac[index,2]
# }

aaa<- sapply(plotData[,1], function(x){
    index<- findInterval(x, ac[,1])
    return(ac[index,2])
})

# aa<- cumsum(tabulate(match(plotData[,1], plotData[,1])))/NROW(plotData)
# aa<- cumsum(sort(plotData[,1]))/sum(plotData[,1])
# plot(aa, sort(plotData[,2]), xlim=c(0,1), ylim=c(0,1)) 

plot(aaa, plotData[,2], xlim=c(0,1), ylim=c(0,1)) 
abline(0,1)
# plot(plotData[,2], ( cumsum(plotData[,1])/sum(plotData[,1]) ),  xlim=c(0,1), ylim=c(0,1)) 
 


findCumProb<- function(freq, cumProb, fetchIndex){
    index<- findInterval(freq, cumProb[,4])
    return( cumProb[index, fetchIndex] )
}


rowSumDataRef<- rowSums(dataRef)
colSumDataRef<- colSums(dataRef)
rowSumDataRefDirty<- rowSums(dataRefDirty)

freqDataRef<- dataRef/rowSumDataRef
freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty




plot(sort(plotData[,2]), aa,  xlim=c(0,1), ylim=c(0,1))
plot(sort(z[1:5000,1]), sort(plotData[,2]),  xlim=c(0,1), ylim=c(0,1)) 

#### old

#simple multinomial (no bias)
p1<- (colSumDataRef[1]+colSumDataRef[2])/2
prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
ll <- sum(log(prob)*colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, paste0("\nMultinomial ", subName2, "\n") )

### end old
























r <- mdmAugmentData(allX)
z <- mdmLogLikeCore(r$y,r$w,f,params)

mdmSingleLogLikeCore(r$y, params[1,])

w<- r$w
r<- r$y
s<- r$y



phi<- 0.01
alpha <- (1-phi)/phi
N <- 10

a <- alpha * 0.2
k <- 3

gamma(alpha)  / gamma(alpha+N)
1/prod( alpha:(alpha+N-1))

# gamma(1-phi) / gamma(1-phi+phi*N) 
1/ prod( (1-phi+phi*(0:(N-1))) ) * phi^N



gamma(alpha)  / gamma(alpha+N) * gamma(a+k) / gamma(a) * gamma(alpha - a + N - k) / gamma(alpha - a)

1 / prod( (1-phi+phi*(0:(N-1))) )  * prod( 0.2*(1-phi)+ phi*(0:(k-1)) ) * prod ( 0.8*(1- phi) + phi*(0:(N-k-1))) 




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
 		norm <- log(sapply(r[,KK], function(x){ prod( (1-phi+phi*(0:(x-1))) ) }))
#  		norm2 <- lgamma( (1-phi) ) - lgamma( (1-phi) +r[,KK])
		pr[,m] <- pr[,m] + norm
	}
	
	ww <- log(rowSums(exp(pr)))
# 	pr <- pr-ww
	list(ll = sum(w*ww), w = exp(pr))
}




# 
# int read_count = data.reads[0]+data.reads[1]+data.reads[2]+data.reads[3];
# 
# 
# 	double alpha_total = alphas[0]+alphas[1]+alphas[2]+alphas[3];
# 	double result = 0.0;
# 	for(int i : {0,1,2,3}) {
# 		for(int x = 0; x < data.reads[i]; ++x) {
# 			result += log(alphas[i]+x);
# 
# 		}
# 	}
# 	for(int x = 0; x < read_count; ++x){
# 		result -= log(alpha_total+x);
# 	}
#
#

condMDM<- function(x, mModel, condIndex=1){

    alpha<- mdmAlphas(mModel$params)
    prop<- mModel$f
    
    nCat<- NROW(alpha)
    
    N<- sum(x)
    
    r<- 0
    
    for(cIndex in 1:nCat){
        for(condIndex in 1:3){
            p<- (N-x[condIndex] + alpha[cIndex, condIndex])/ (N + sum(alpha[cIndex, ]) -1  )
            cat(p, "\t")
            r <- r+p
        }
        cat(r, "\n")
        r<- 0
    }

    
    
    
}
