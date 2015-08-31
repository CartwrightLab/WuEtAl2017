#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

##########################
###
##########################

####################
#### from don
# # 
# # Hi Rachel,
# # 
# # I'm attaching the two lists that you sent me annotated with 4
# # different overlaps. each overlap is binary (0=no overlap, 1=overlap).
# # I assumed your coordinates are w.r.t. to hg19/NCBI37
# # 
# # -sds = overlaps a sequence with a large (>1kb), highly homologous
# # (>90% identity) duplication somewhere else in the genome
# # -cnv=overlaps a CNV from Conrad, et al 2009
# # -cent=within 1Mb of a centromere
# # -tel=within 1Mb of a telomere
# # 
# # Looks like you have a lot of hits for sds and cnvs.
# # 
# # Cheers
# # 
# # Don
##############################

# pid <- Sys.getpid()
# args <- commandArgs(trailingOnly=TRUE)
dirty_data <- FALSE

# nn <- as.numeric(args[1])
# k <- args[2]
# dirty_data <- grepl("^\\d+[Dd]$",k)
# if(dirty_data) {
# 	k <- as.numeric(sub("[Dd]$","",k))
# } else {

pwd <- "/home/steven/Postdoc2/Project_MDM/run/"
setwd(pwd)

load("summaryResult2011.RData")
sites878<- read.table("/home/steven/Postdoc2/Project_MDM/run/sites878_labeled.txt", header=T, comment.char="")

index<- sites878[,2] %in% rownames(dataRef)
sitesSub<- sites878[index,]
if(!all.equal(as.character(sitesSub[,2]) , rownames(dataRef))){
	stop("Different/mismatched names: sitesSub, rownames(dataRef)")
}

index<- sites878[,2] %in% rownames(dataFP)
sitesSubFP<- sites878[index,]
if(!all.equal(as.character(sitesSubFP[,2]) , rownames(dataFP))){
	stop("Different/mismatched names: sitesSub, rownames(dataFP)")
}



apply(maxLikelihoodTable[[3]],1,which.max)


index<- apply(maxLikelihoodFP[[3]],1,which.max)
table(sitesSubFP[index==1,3:4])
table(sitesSubFP[index==2,3:4])

###########################################
##### from summary.R
subFolders <- c(
# "2011/multiallelic/base_count/", 
"2011/original/base_count/"
)
p<- 1
	
	fullPath <- paste(pwd, subFolders[[p]], sep="")


	setwd(fullPath)
	maxModel<- extractMaxModel(fullPath)
	

	upperLimit <- 110
	lowerLimit <- 20
	dirtyData <- F
	
	hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
	dat <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
	dataRef<- parseData(dat, lowerLimit, upperLimit, dirtyData)
	
	sapply(maxModel, function(x){x$ll})
	sapply(maxModel, function(x){x$f})
	sapply(maxModel, function(x){x$params})

	maxLikelihoodTable<- calculateEachLikelihood(dataRef, maxModel)#, numData=100)

	modelParametersSummary<- sapply(maxModel, function(x){
		x<- cbind(x$params, mdmAlphas(x$params), prop=x[["f"]] )
		x<- x[,c(2:7,1,8)]
	})

propEM<-sapply(maxModel, function(x){x$f})
propML<- getMaxLikelihoodProp(maxLikelihoodTable)

modelParametersSummary2 <- modelParametersSummary
propResultTable<- vector("list", length(propEM))
names(propResultTable)<- names(propEM)
for(i in 1:length(propML)){
	diff <- propEM[[i]] - propML[[i]]
	propResultTable[[i]]<- rbind(propEM[[i]], propML[[i]], diff)
	rownames(propResultTable[[i]])<- c("EM", "ML", "Diff")
	
	if(i < 3){
		modelParametersSummary2[[i]] <- cbind(t(modelParametersSummary[[i]]), "ML"=propML[[i]], "Diff"=diff)
	}
	else{
		modelParametersSummary2[[i]] <- cbind(modelParametersSummary[[i]], "ML"=propML[[i]], "Diff"=diff)
	}
}


totalAbsDiff<- vector()
for(i in 1:length(propML)){
	totalAbsDiff[i]<- sum(abs(propML[[i]] - propEM[[i]]))/2
	
}




summary(apply(maxLikelihoodTable[[2]], 1, diff))
plot(sort(apply(maxLikelihoodTable[[2]], 1, diff) ))

apply(maxLikelihoodTable[[2]],1,which.max)

prop/sum(prop)
maxModel[[d]]$f 

propML
propEM

params<- maxModel[[1]]$params




###################
## check files

allFiles<- list.files(pattern="*byref*")

dat<- list()
count <- list()
for(f in 1:length(allFiles)){
	dat[[f]] <- read.delim(allFiles[f],header=TRUE)
	count[[f]] <- (dat[[f]]$snp == 1 & dat[[f]]$snpdif == 0)
# 	x <- cbind(dat$refs,dat$alts,dat$e1s+dat$e2s)
# 	x <- data.matrix(x)
# 	row.names(x) <- dat$pos
# 	x <- x[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirty_data), ]
# 	n <- rowSums(x)
# 	oo <- lowerLimit <= n & n <= upperLimit

}
# x <- x[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirty_data), ]


> getwd()
[1] "/home/steven/Postdoc2/Project_MDM/run/2010/original/base_count"
> sapply(count,sum)
[1] 28979 18645     0
> sapply(dat,dim)
      [,1]  [,2]     [,3]
[1,] 37009 20207 35023641
[2,]    11    11       11









