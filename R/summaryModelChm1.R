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

# pid <- Sys.getpid()
# args <- commandArgs(trailingOnly=TRUE)

# nn <- as.numeric(args[1])
# k <- args[2]
# dirtyData <- grepl("^\\d+[Dd]$",k)
# if(dirtyData) {
# 	k <- as.numeric(sub("[Dd]$","",k))
# } else {
#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

# source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/mdm.R")

isCEU <- FALSE


dirtyData <- TRUE
upperLimit <- 150
lowerLimit <- 10
loadData<- TRUE


latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"


if(!isCEU){
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
}
setwd(pwd)


p<- 2
p<- 1

# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

setwd(fullPath)
maxModel<- extractMaxModel(fullPath, isCEU)

if(!isCEU){
    hets_byref<- file.path("base_count_meta_subsample") 
}


dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU=isCEU)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU=isCEU)
n <- rowSums(dataRefDirty)
propRef <- dataRefDirty[,1]/n
oo <- propRef > 0.8
dataRefProp <- dataRefDirty[oo,]


fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTableFull.RData")
if ( file.exists(fileMaxLikelihoodTabel) && loadData ){
    load(fileMaxLikelihoodTabel)
} else{
    maxLikelihoodTable<- calculateEachLikelihoodCHM1(maxModel, dataFull, lowerLimit=lowerLimit, upperLimit=upperLimit, isCEU=isCEU)#, numData=100)
    attr(maxLikelihoodTable, "title")<- subName
    save(maxLikelihoodTable, file=file.path(fullPath, "maxLikelihoodTableFull.RData"))
}

modelLikelihood<- sapply(maxModel, function(x){x$ll})
# modelParams<- sapply(maxModel, function(x){x$params})

modelParametersSummary<- sapply(maxModel, function(x){
    x<- cbind(x$params, mdmAlphas(x$params), prop=x[["f"]] )
    x<- x[,c(2,4,5,7,1,8)]
    x<- matrix(x, ncol=6)
    colnames(x)<- c("pi_ref", "pi_err", "alpha_ref", "alpha_err", "phi", "prop")
    
    return(x)
})
modelParametersSummary<- modelParametersSummary[grepl("_[0-9]D|P",names(modelParametersSummary))]



propEM<- sapply(maxModel, function(x){ x$f })
propEM<- propEM[grepl("_[0-9]D|P",names(propEM))]
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


#     newOrder<-rev(order(modelParametersSummary2[[i]][,8]))
#     modelParametersSummary2[[i]]<- 
#         modelParametersSummary2[[i]][newOrder,]
    
}


modelParametersSummary3<- vector(length=length(modelParametersSummary2)+1, mode="list")
modelParametersSummary3[2:13] <- modelParametersSummary2
names(modelParametersSummary3)<- c("Multinomial", names(modelParametersSummary2))

colSumDataRef<- colSums(dataRefProp)[-2]
#simple multinomial (no bias)
prob_MN <- colSumDataRef /sum(colSumDataRef)
modelParametersSummary3[[1]] <- matrix(c(prob_MN,rep(NA,5)), nrow=1, byrow=T)
# names(modelParametersSummary3[[1]])<- colnames(modelParametersSummary3[[3]])


## latex parameter table    
# formatC(x, digits=3, width=8, format="g")
# sprintf("%8.3g" , x)
latexTable<- sapply(modelParametersSummary3, function(x){
    newOrder<-rev(order(x[,6]))
    x<- x[newOrder,]
    
    x<- formatC(x, digits=3, width=8, format="g")
    x<- matrix(x,ncol=7)
    s<- apply(x,1,function(y){
        paste(y, collapse=" & ")
    })
    latex<- paste(s, collapse=" \\\\ \n & ")
})



header<- names(latexTable)
# header<- gsub("hets_" , "", header)
# header<- gsub("_" , "M", header)
# 
# header<- paste(subName, 1:length(latexTable), sep=" M" )
header<- paste("CHM1", header)
header<- gsub("_" , "M", header)
header<- gsub("D" , "F", header)
header<- gsub("P" , "R", header) 
 
prefix<- paste0("\\begin{tabular}{|c|cc|cc|c|c|c|}
    \\hline \\multicolumn{8}{|c|}{Parameter estimates ",fullTitle,"}\\\\ \\hline
    Model & $\\pi_{ref}$ & $\\pi_{err}$ & $\\alpha_{ref}$ & $\\alpha_{err}$ & $\\varphi$ &  $p$ & ML-P \\\\ \\hline")

sufix<- "\\end{tabular}"

latexTableFull <- latexTable
fileLatexTable <- file.path(latexDir, paste0(subName,"_latexTable.tex") )
cat(prefix, file=fileLatexTable, fill=T)
for(i in 1:length(latexTable) ){
    latexTableFull[i]<- paste0(header[i], " & ", latexTable[[i]], " \\\\ \\hline")
    cat(latexTableFull[i], file=fileLatexTable, fill=T, append=T)
}
cat(sufix, file=fileLatexTable, fill=T, append=T)
# write.table(latexTableFull, file=file.path(fullPath, "latxTable"), quote=F, row.names=F)
    

## maxlikelihood table
modelLikelihood<- sapply(maxModel, function(x){x$ll})
modelLikelihood<- modelLikelihood[grepl("_[0-9]D|P",names(modelLikelihood))]

# \begin{tabular}{|c|c|c|}
#         \hline \multicolumn{3}{|c|}{-lnL}\\ \hline
#         Model & filtered & unfiltered \\\hline
#         DM & -1457308 & -1623316 \\ 
#         2 DM & -1456421 &-1619317  \\ 
#         3 DM & -1456400 & -1619209   \\ 
#         4 DM &-1456380  &  -1619184  \\
#         5 DM &-1456372 & -1619176  \\
#         6 DM & -1456368& -1619169  \\\hline
#       \end{tabular}


fileMaxLikelihoodLatexTabel <- file.path(latexDir, paste0(subName, "_maxLikelihoodLatexTable.tex") )

prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|c|c|}
    \\hline \\multicolumn{7}{|c|}{",fullTitle," } \\\\ \\hline
     Model & lnL & lnL P & AIC & BIC & AIC P& BIC P\\\\ \\hline")
sufix<- "\\hline\n\\end{tabular}"
#      Model & lnL & AIC & BIC \\\\ \\hline")
# maxLiTable<- matrix(ncol=3, nrow=length(modelLikelihood) )
# numFreeP<- seq(3,by=4,length=6)
# coefBIC<- log(NROW(dataRefDirty)) 
# for(i in 1:NROW(maxLiTable) ){
#     t1<- modelLikelihood[i]
#     t2<- -2*t1 + 2*numFreeP[i]
#     t3<- -2*t1 + coefBIC*numFreeP[i]
#     maxLiTable[i,]<- c(t1, t2, t3)
# }
maxLiTable<- matrix(ncol=6, nrow=length(modelLikelihood)/2 )
numFreeP<- seq(3,by=4,length=6)
coefBIC<- c(log(NROW(dataRefDirty)), log(NROW(dataRefProp)) )
for(i in 1:NROW(maxLiTable) ){
    t1<- c(modelLikelihood[i*2-1], modelLikelihood[i*2])
    t2<- -2*t1 + 2*numFreeP[i]
    t3<- -2*t1 + coefBIC*numFreeP[i]
    maxLiTable[i,]<- c(t1, t2, t3)
}


minIndex<- apply(maxLiTable,2,which.min)
maxLiTable<- formatC(maxLiTable,  digits=2, width=12, format="f")
for(i in 3:NCOL(maxLiTable)){
    maxLiTable[minIndex[i], i]<-paste0(maxLiTable[minIndex[i], i], "*")
}

maxLikelihoodLatex<- apply(maxLiTable,1,function(y){
        paste(y, collapse=" & ")
    })


cat(prefix, file=fileMaxLikelihoodLatexTabel, fill=T)
for(i in 1:length(maxLikelihoodLatex) ){
    latex<- paste0(header[i], " & ", maxLikelihoodLatex[i], " \\\\ ")
    cat(latex, file=fileMaxLikelihoodLatexTabel, fill=T, append=T)
}
cat(sufix, file=fileMaxLikelihoodLatexTabel, fill=T, append=T)




# }  ## end for loop


    
#writeLines 
# totalAbsDiff<- vector()
# for(i in 1:length(propML)){
#     totalAbsDiff[i]<- sum(abs(propML[[i]] - propEM[[i]]))/2
# 
# }
#####################################################################################
##### NOT YET UPDATE AFTER THIS POINT
############################################################
#######################################
##### FP. 
###############################
# dataRef<- parseData(dat, lowerLimit, upperLimit, dirtyData)


p<- 8
# for(p in 1:length(subNameList) ){

# subName<- subNameList[p]
# fullTitle<- fullTitleList[p]
# subFolders <- paste0(subName, "/original/base_count/")
# fullPath <- file.path(pwd, subFolders)
# 
# setwd(fullPath)
# maxModel<- extractMaxModel(fullPath)
# 
# hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
# dat <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dat, lowerLimit, upperLimit, dirtyData)
# dataRefDirty<- parseData(dat, lowerLimit, upperLimit, dirtyData=TRUE)

## SNP count
if(isCEU){
    fileSnpCountLatexTabel <- file.path(latexDir, paste0("snpCountLatexTable.tex") )
} else {
    fileSnpCountLatexTabel <- file.path(latexDir, paste0("snpCountLatexTable_CHM1.tex") )
}


prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|}
    \\hline \\multicolumn{5}{|c|}{", "" ," } \\\\ \\hline
    Dataset & Method (1) & Method (2) & Both methods & True heterozygotes \\\\ \\hline")
sufix<- "\\hline\n\\end{tabular}"


cat(prefix, file=fileSnpCountLatexTabel, fill=T)
for(p in length(subNameList):1 ){

    subName<- subNameList[p]
    fullTitle<- fullTitleList[p]
    subFolders <- paste0(subName, "/original/base_count/")
    fullPath <- file.path(pwd, subFolders)
    setwd(fullPath)

    hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
    dat <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)

    name2<- gsub("_" , " ", subName)
    
    if(isCEU){
        snpCount <- c(name2, table(dat$callby)[c(1,3,2)], sum( (dat$callby == 2 & (dat$snp == 1 & dat$snpdif == 0)) ) )  
    } else {
        snpCount <- c(name2, table(dat$callby)[1], 0,0, sum( (dat$callby == 1 & (dat$snp == 1 & dat$snpdif == 0)) ) )  
    }
 
    snpLatex <- paste(snpCount, collapse=" & ")
    snpLatex <- paste(snpLatex ,"\\\\")
    cat(snpLatex, file=fileSnpCountLatexTabel, fill=T, append=T)
}
cat(sufix, file=fileSnpCountLatexTabel, fill=T, append=T)



####################



x <- x[dat$callby == 2 & (dat$snp == 1 & dat$snpdif == 0) , ]

indexCall2Snp <- (dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 0) ))
indexCall2SnpDiff <- (dat$callby == 2 & ((dat$snp == 1 & dat$snpdif == 1) ))


table(dat$callby, dat$snp, dat$snpdif)

table(dat$callby)


table(dat$callby, dat$snp)

table(dat$callby, dat$snpdif)

table(dat$snp, dat$snpdif)

indexFP<- dat$callby == 2 & dat$snp == 0
indexTP <- dat$callby == 2 & dat$snp == 1


dataFP<- parseDataIndex(dat, indexFP, lowerLimit, upperLimit)
dataTP<- parseDataIndex(dat, indexTP, lowerLimit, upperLimit)
dataSNP<- parseDataIndex(dat, indexCall2Snp, lowerLimit, upperLimit)
dataSNP<- parseDataIndex(dat, indexCall2SnpDiff, lowerLimit, upperLimit)

dataXP<- list(FP=dataFP, TP=dataTP)
pp<- list()
for(x in 1:length(dataXP)){
    pp[[x]]<- list()
    for(d in 1:length(maxModel)){
        maxLikelihoodXP<- calculateEachLikelihoodOneModel(maxModel[[d]], dataXP[[x]])
        pp[[x]][[d]]<- apply(maxLikelihoodXP,1,which.max)
    }
    pp[[x]]
}
names(pp)<- names(dataXP)

maxModelIndex<- sapply(maxModel,function(x){which.max(x$f)})
maxModelIndex
sapply(pp$FP,function(x){which.max(table(x))})
sapply(pp$TP,function(x){which.max(table(x))})


sapply(pp$FP,function(x){
    y<- (table(x))
    s<- sum(y)
    r<- y/s
    rownames(r)<- NULL
    return(r)
})


sapply(pp$TP,function(x){
    y<- (table(x))
    s<- sum(y)
    r<- y/s
    rownames(r)<- NULL
    return(r)
})


# > p
# [1] 4
# > maxModelIndex
#     hets_CEU13_1    hets_CEU13_1D  hets_CEU13_2.m2 hets_CEU13_2D.m2 
#                1                1                2                2 
#  hets_CEU13_3.m3 hets_CEU13_3D.m3  hets_CEU13_4.m4 hets_CEU13_4D.m4 
#                3                3                4                4 
#  hets_CEU13_5.m5 hets_CEU13_5D.m5  hets_CEU13_6.m6 hets_CEU13_6D.m6 
#                5                5                6                6 
# > sapply(pp$FP,function(x){which.max(table(x))})
# 1 1 1 1 1 2 1 2 2 4 5 5 
# 1 1 1 1 1 2 1 2 2 4 5 5 
# > sapply(pp$TP,function(x){which.max(table(x))})
# 1 1 2 2 3 3 4 4 5 5 6 6 
# 1 1 2 2 3 3 4 4 5 5 6 6 
# 

# > sapply(pp$FP,function(x){
# +     y<- (table(x))
# +     s<- sum(y)
# +     r<- y/s
# +     rownames(r)<- NULL
# +     return(r)
# + })
# [[1]]
# [1] 1
# 
# [[2]]
# [1] 1
# 
# [[3]]
# [1] 0.6591884 0.3408116
# 
# [[4]]
# [1] 0.6726053 0.3273947
# 
# [[5]]
# [1] 0.66336743 0.02155504 0.31507753
# 
# [[6]]
# [1] 0.2096118 0.4480370 0.3423513
# 
# [[7]]
# [1] 0.47740020 0.12537116 0.03156274 0.36566590
# 
# [[8]]
# [1] 0.28571429 0.37611349 0.02397449 0.31419773
# 
# [[9]]
# [1] 0.03046299 0.45353569 0.11613329 0.05454745 0.34532058
# 
# [[10]]
# [1] 0.17211041 0.27361707 0.03079292 0.28978335 0.23369625
# 
# [[11]]
# [1] 0.08567030 0.03145277 0.05036842 0.12669086 0.46376333 0.24205433
# 
# [[12]]
# [1] 0.26954800 0.05971627 0.03024304 0.13691851 0.29385241 0.20972176
# 
# > 
# > 
# > sapply(pp$TP,function(x){
# +     y<- (table(x))
# +     s<- sum(y)
# +     r<- y/s
# +     rownames(r)<- NULL
# +     return(r)
# + })
# [[1]]
# [1] 1
# 
# [[2]]
# [1] 1
# 
# [[3]]
# [1] 0.185236 0.814764
# 
# [[4]]
# [1] 0.1850017 0.8149983
# 
# [[5]]
# [1] 0.212956143 0.004117844 0.782926013
# 
# [[6]]
# [1] 0.09782390 0.09112822 0.81104787
# 
# [[7]]
# [1] 0.117609642 0.073284232 0.005390023 0.803716103
# 
# [[8]]
# [1] 0.154871108 0.051054570 0.006494811 0.787579511
# 
# [[9]]
# [1] 0.002477402 0.097556076 0.150284566 0.011416137 0.738265819
# 
# [[10]]
# [1] 0.15477067 0.20984265 0.01158353 0.01797790 0.60582524
# 
# [[11]]
# [1] 0.083126883 0.002912621 0.010010044 0.304854369 0.105724807 0.493371276
# 
# [[12]]
# [1] 0.20954135 0.01382658 0.01161701 0.21118179 0.01817877 0.53565450
# 
# > 


table(ppFP)
# pp
#    1    2 
# 5994 3099 
maxLikelihoodTP<- calculateEachLikelihoodOneModel(maxModel[[3]], dataTP)
ppTP<- apply(maxLikelihoodTP,1,which.max)
table(ppTP)


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









