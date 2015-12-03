#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

isCEU <- FALSE
# isCEU <- TRUE


dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10

latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

# if(isCEU){
#     pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
#     subNameList<- c(
#     "CEU10_C10", "CEU10_C21",
#     "CEU11_C10", "CEU11_C21",
#     "CEU12_C10", "CEU12_C21",
#     "CEU13_C10", "CEU13_C21"
#     )
# 
#     fullTitleList<- c(
#     "CEU 2010 Chromosome 10", "CEU 2010 Chromosome 21",
#     "CEU 2011 Chromosome 10", "CEU 2011 Chromosome 21",
#     "CEU 2012 Chromosome 10", "CEU 2012 Chromosome 21",
#     "CEU 2013 Chromosome 10", "CEU 2013 Chromosome 21"
#     )
# } else {
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
# }
setwd(pwd)


collapseSortMean<- function(data, ncol){
    result<- apply(data, 2, function(x){
        xm<- matrix(x, ncol=ncol)
        xmSort<- apply(xm, 2, function(y){ sort(y) } )
        xmMean<- apply(xmSort, 1, function(y){ mean(y) } )
        return(xmMean)

    })
    return(result)
}


plotNamePrefix<- "qqPlotsV2_"
numRep<- 100


p<- 2
p<- 1

# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# # setwd(fullPath)
# hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
# dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU)
# dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

hets_byref<- file.path("base_count_meta_subsample") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU=isCEU)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU=isCEU)
n <- rowSums(dataRefDirty)
propRef <- dataRefDirty[,1]/n
oo <- propRef > 0.8
# oo <- lowerLimit <= n & n <= upperLimit 
dataRefProp <- dataRefDirty[oo,]

plotData<- dataRefDirty[,-2]
plotData<- dataRefProp[,-2]
rowSumDataRef<- rowSums(plotData)
colSumDataRef<- colSums(plotData)
# rowSumDataRefDirty<- rowSums(dataRefDirty)

freqDataRef<- plotData/rowSumDataRef
# freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty
subName2<- gsub("_C", " Chr", subName)

maxModel<- extractMaxModel(fullPath, isCEU=isCEU)
whichIsDirty <- grepl("_[0-9]D",names(maxModel))
whichIsP <- grepl("_[0-9]P",names(maxModel))
header<- gsub("_" , "", names(maxModel) )
header<- gsub("D" , " FD", header )
header<- gsub("P" , " components", header ) #need this one
header<-paste(subName2, header)

plotTitle <- paste0(plotNamePrefix, subName, ".pdf") 
qqplotFile<- file.path(latexDir, plotTitle)


pdf(file=qqplotFile, width=12, height=6, title=plotTitle)
par(mai=c(0.8,0.8,0.2,0.1), mfrow=c(1,2), 
    cex.main=1.2^2,cex.lab=1.2, 
    omi=c(0,0,0.5,0) )


## simple multinomial (no bias)
prob <- colSumDataRef /sum(colSumDataRef)
ll <- sum(log(prob)*colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
b2 <- b/rowSums(b)
expFreq<- collapseSortMean(b2, numRep)

plotqq(expFreq, freqDataRef, paste0("\nMultinomial ", subName2, "\n") )


 
# ## multinomial (with ref bias)
# prob <- colSumDataRef / sum(colSumDataRef)
# ll <- sum(log(prob)* colSumDataRef)
# # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
# b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
# b2 <- b/rowSums(b)
# expFreq<- collapseSortMean(b2, numRep)
# plotqq(expFreq, freqDataRef, paste0("\nBiased Multinomial ", subName2, " \n"))

## dirichlet-multinomial mixture
for( m in which(whichIsP) ) {

    main<- paste0("\nMixture of Dirichlet Multinomial ", header[m] ," \n")
#     cat("Plotting ", main)
    mModel<- maxModel[[m]]
    
    b <- rmdm(numRep*length(rowSumDataRef), rep(rowSumDataRef, numRep), mModel$f,
                    phi=mModel$params[,1], p=mModel$params[,c(2,4)])
    b2 <- b/rowSums(b)
    expFreq<- collapseSortMean(b2, numRep)
    
    main<- paste0("\nMixture of Dirichlet Multinomial ", header[m] ," \n")
    if(m==3){
        main<- gsub(" 1 components", "", main)
        main<- gsub("Mixture of ", "", main)
    }
    
    plotqq(expFreq, freqDataRef, main)
}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")
# } # match (p in 1:length(subNameList) ){

