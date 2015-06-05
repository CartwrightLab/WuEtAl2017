##########################
### TODO:
### [x] Rerun EM for 1000 iterations
### [x] Try a different chromosome, or a section
### Re-estimate categories after divided into good/bad sites.
###		Hopefully bad sites are at minor categories
### 
### Almost done, just need to clean up. Get estimate for all sites (from 4 datasets) 
###
##########################

#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10


# pid <- Sys.getpid()
# args <- commandArgs(trailingOnly=TRUE)

# nn <- as.numeric(args[1])
# k <- args[2]
# dirtyData <- grepl("^\\d+[Dd]$",k)
# if(dirtyData) {
# 	k <- as.numeric(sub("[Dd]$","",k))
# } else {

pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
setwd(pwd)

subNameList<- c(
"CEU12_C10",
"CEU12_C21",
"CEU13_C10",
"CEU13_C21"
)



p<- 4
# for(p in 1:length(subName) ){

subName<- subNameList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

setwd(fullPath)
maxModel<- extractMaxModel(fullPath)

hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dat <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
dataRef<- parseData(dat, lowerLimit, upperLimit, dirtyData)

fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTable.RData")
if ( file.exists(fileMaxLikelihoodTabel) ){
    load(fileMaxLikelihoodTabel)
} else{
#     maxLikelihoodTable<- calculateEachLikelihood(dataRef, maxModel)#, numData=100)
    attr(maxLikelihoodTable, "title")<- subName
    save(maxLikelihoodTable, file=file.path(fullPath, "maxLikelihoodTable.RData"))
}

modelLikelihood<- sapply(maxModel, function(x){x$ll})
# modelParams<- sapply(maxModel, function(x){x$params})

modelParametersSummary<- sapply(maxModel, function(x){
    x<- cbind(x$params, mdmAlphas(x$params), prop=x[["f"]] )
    x<- x[,c(2:7,1,8)]
    x<- matrix(x, ncol=8)
    colnames(x)<- c("pi_ref", "pi_alt", "pi_err", "alpha_ref", "alpha_alt", "alpha_err", "phi", "prop")
    return(x)
})




propEM<- sapply(maxModel, function(x){x$f})
propML<- getMaxLikelihoodProp(maxLikelihoodTable)

modelParametersSummary2 <- modelParametersSummary
propResultTable<- vector("list", length(propEM))
names(propResultTable)<- names(propEM)
for(i in 1:length(propML)){
    diff <- propEM[[i]] - propML[[i]]
    propResultTable[[i]]<- cbind(propEM[[i]], propML[[i]], diff)
    colnames(propResultTable[[i]])<- c("EM", "ML", "Diff")
    
#     if(i < 3){
#         modelParametersSummary2[[i]] <- cbind(t(modelParametersSummary[[i]]), "ML"=propML[[i]], "Diff"=diff)
#     } else{
        modelParametersSummary2[[i]] <- cbind(modelParametersSummary[[i]], "ML"=propML[[i]])
#         "Diff"=diff)
#     }
}
    
# formatC(x, digits=3, width=8, format="g")
# sprintf("%8.3g" , x)
latexTable<- sapply(modelParametersSummary2, function(x){
    x <- formatC(x, digits=3, width=8, format="g")
    s<- apply(x,1,function(y){
        paste(y, collapse=" & ")
    })
    latex<- paste(s, collapse=" \\\\ \n & ")
})

header<- names(latexTable)
header<- gsub("hets_" , "", header)
header<- gsub("_" , "M", header)
latexTableFull <- latexTable
fileLatexTable <- file.path(pwd, paste0(subName,"_latexTable.tex") )
cat("",file=fileLatexTable)
for(i in 1:length(latexTable) ){
    latexTableFull[i]<- paste0(header[i], " & ", latexTable[[i]], " \\\\ \\hline")
    cat(latexTableFull[i],file=fileLatexTable, fill=T, append=T)
}
# write.table(latexTableFull, file=file.path(fullPath, "latxTable"), quote=F, row.names=F)
    

#writeLines

# totalAbsDiff<- vector()
# for(i in 1:length(propML)){
#     totalAbsDiff[i]<- sum(abs(propML[[i]] - propEM[[i]]))/2
# 
# }

#######################################
##### FP. in child not NOT in human
###############################
x <- x[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]

table(dat$callby, dat$snp, dat$snpdif)

table(dat$callby)


data.fp <- dat[dat$callby == 2 & dat$snp == 0, ]


dataFP<- parseDataFP(dat, lowerLimit, upperLimit)

maxLikelihoodFP<- calculateEachLikelihood(dataFP, maxModel)

propFP<- getMaxLikelihoodProp(maxLikelihoodFP)


cbind(unlist(propEM), unlist(propML), unlist(propFP) )

modelParametersSummary2 <- modelParametersSummary
propResultTable<- vector("list", length(propEM))
names(propResultTable)<- names(propEM)
for(i in 1:length(propML)){
	diff_ML <- propEM[[i]] - propML[[i]]
	diff_FP <- propEM[[i]] - propFP[[i]]
	propResultTable[[i]]<- rbind(propEM[[i]], propML[[i]], propFP[[i]], diff_ML, diff_FP)
	rownames(propResultTable[[i]])<- c("EM", "ML", "FP", "Diff_ML", "Diff_FP")
	
	if(i < 3){
		modelParametersSummary2[[i]] <- cbind(t(modelParametersSummary[[i]]), "ML"=propML[[i]], "FP"=propFP[[i]], "Diff_ML"=diff_ML, "Diff_FP"=diff_FP)
	}
	else{
		modelParametersSummary2[[i]] <- cbind(modelParametersSummary[[i]], "ML"=propML[[i]], "FP"=propFP[[i]], "Diff_ML"=diff_ML, "Diff_FP"=diff_FP)
	}
}

####################


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
# 	x <- x[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]
# 	n <- rowSums(x)
# 	oo <- lowerLimit <= n & n <= upperLimit

}
# x <- x[dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) | dirtyData), ]


> getwd()
[1] "/home/steven/Postdoc2/Project_MDM/run/2010/original/base_count"
> sapply(count,sum)
[1] 28979 18645     0
> sapply(dat,dim)
      [,1]  [,2]     [,3]
[1,] 37009 20207 35023641
[2,]    11    11       11









