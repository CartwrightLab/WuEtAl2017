#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

isCEU <- FALSE
isCEU <- TRUE

loadData <- TRUE
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


plotNamePrefix<- "qqPlotsV2_"
numRep<- 100

p<- 8

for(p in 1:length(subNameList) ){

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

qqplotFile<- file.path(latexDir, paste0(plotNamePrefix, subName, ".pdf") )


pdf(file=qqplotFile, width=15, height=6, title=qqplotFile)
par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), 
    cex.main=1.2^4,cex.lab=1.2^2, 
    omi=c(0,0,0.5,0) )

# mains = c("Reference Allele", "Alternate Allele", "Error")

## simple multinomial (no bias)
plotName<- "Multinomial"
fileExpectedFreq <- paste0(getwd(), "/", subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData")

if ( file.exists(fileExpectedFreq) && loadData ){
    load(fileExpectedFreq)
} else{
    p1<- (colSumDataRef[1]+colSumDataRef[2])/2
    prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
    ll <- sum(log(prob)*colSumDataRef)
    # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
    b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
    b2 <- b/rowSums(b)
    expFreq<- collapseSortMean(b2, numRep)
    save(expFreq, file=fileExpectedFreq)
}
plotqq(expFreq, freqDataRef, paste0(plotName, " ", subName2) )


## multinomial (with ref bias)
plotName<- "Biased Multinomial"
fileExpectedFreq <- paste0(getwd(), "/", subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData")

if ( file.exists(fileExpectedFreq) && loadData ){
    load(fileExpectedFreq)
} else{

    prob <- colSumDataRef / sum(colSumDataRef)
    ll <- sum(log(prob)* colSumDataRef)
    # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
    b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
    b2 <- b/rowSums(b)
    expFreq<- collapseSortMean(b2, numRep)
    save(expFreq, file=fileExpectedFreq)
}
plotqq(expFreq, freqDataRef, paste0(plotName, " ", subName2) )


## dirichlet-multinomial mixture
for( m in which(! whichIsDirty) ) {
    
    plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
    if(m==1){
        plotName<- gsub(" 1 components", "", plotName)
        plotName<- gsub("Mixture of ", "", plotName)
    }

    fileExpectedFreq <- paste0(getwd(), "/", subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData")

    if ( file.exists(fileExpectedFreq) && loadData ){
        load(fileExpectedFreq)
    } else{
        mModel<- maxModel[[m]]
        b <- rmdm(numRep*length(rowSumDataRef), rep(rowSumDataRef, numRep), mModel$f,
                        phi=mModel$params[,1], p=mModel$params[,2:4])
        b2 <- b/rowSums(b)
        expFreq<- collapseSortMean(b2, numRep)
    
        save(expFreq, file=fileExpectedFreq)
    }
    plotqq(expFreq, freqDataRef, plotName)
}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")
} # match (p in 1:length(subNameList) ){


####################################################################
##### Plot MS
####################################################################

qqplotFile<- file.path(latexDir, paste0(plotNamePrefix, "CEU_MS.pdf") )
pdf(file=qqplotFile, width=15, height=6, title=qqplotFile)
par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), 
    cex.main=1.2^4,cex.lab=1.2^2, 
    omi=c(0,0,0,0) )

## simple multinomial (no bias)
plotName<- "Multinomial"
fileExpectedFreq <- paste0(getwd(), "/", subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData")

if ( file.exists(fileExpectedFreq) && loadData ){
    load(fileExpectedFreq)
} else{
    warning(paste("File doesn't exist:", fileExpectedFreq))
}
plotqq(expFreq, freqDataRef, "")


## dirichlet-multinomial mixture
for( m in which(! whichIsDirty)[c(1,3)] ) {
    
    plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
    if(m==1){
        plotName<- gsub(" 1 components", "", plotName)
        plotName<- gsub("Mixture of ", "", plotName)
    }

    fileExpectedFreq <- paste0(getwd(), "/", subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData")

    if ( file.exists(fileExpectedFreq) && loadData ){
        load(fileExpectedFreq)
    } else{
        warning(paste("File doesn't exist:", fileExpectedFreq))
    }
    plotqq(expFreq, freqDataRef, "")
}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")

