#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- FALSE
suppressPackageStartupMessages( source("./summaryFunctions.R") )
source("./summarySetup.R")

plotLabels = function(title) {
    op <- par("mai")
    par(mai=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,title,cex=1.2^4,font=2,xpd=NA)
    par(mai=c(op[1],0,op[3],0.1))
    plot.new()
    text(0.5,0.5,"Observed Frequency",srt=90,cex=1.2^3,xpd=NA)
    par(mai=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,"Expected Frequency",cex=1.2^3,xpd=NA)
    par(mai=op)
}


plotNamePrefix<- "qqPlotsExp_"
numRep<- 100


for(p in length(subNameList):1 ){

    subName<- subNameList[p]
    fullTitle<- fullTitleList[p]
    fullPath <- file.path(dataDir, subName)
   
    temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
    dataFull <- temp$dataFull
#     dataRef <- temp$dataRef
    dataRefDirty <- temp$dataRefDirty

    maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)

    n <- rowSums(dataRefDirty)
    propRef <- dataRefDirty[,1]/n
    oo <- propRef > 0.8
    dataRefProp <- dataRefDirty[oo,]

    plotData<- dataRefDirty[,-2]
    plotData<- dataRefProp[,-2]
    rowSumDataRef<- rowSums(plotData)
    colSumDataRef<- colSums(plotData)

    freqDataRef<- plotData/rowSumDataRef
    subName2<- gsub("_C", " Chr", subName)


    whichIsP <- grepl("_[0-9]P",names(maxModel))
    header<- gsub("_" , "", names(maxModel) )
    header<- gsub("D" , " FD", header )
    header<- gsub("P" , " components", header ) #need this one
    header<-paste(subName2, header)

    plotTitle <- paste0(plotNamePrefix, subName, ".pdf") 
    qqplotFile<- file.path(figureDir, plotTitle)

    mtop = 0.4
    mbottom = 0.3
    mleft = 0.4
    mright = 0
    width = 7.5
    sz = (width-mleft-mright)/3

    # rescale to 1x2 
    height = mtop+mbottom+sz
    width = 2*sz+mleft+mright

    m = c(
        0,1,1,0,
        2,4,5,0,
        0,3,3,0
    );
    m = matrix(m,3,4,byrow=TRUE)
    w = c(mleft,sz,sz,mright)
    h = c(mtop,sz,mbottom)

    cat("Creating File ", qqplotFile, "...\n")
    pdf(file=qqplotFile, width=width, height=height, title=subName, useDingbats=T,pointsize=12)
    par(mai=c(0.3,0.15,0.0,0.15))
    layout(m,w,h)

## simple multinomial (no bias)
plotName<- "Multinomial"
fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )

if ( file.exists(fileExpectedFreq) && loadData ){
    load(fileExpectedFreq)
} else{
    prob <- colSumDataRef /sum(colSumDataRef)
    ll <- sum(log(prob)*colSumDataRef)
    # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

    b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
    b2 <- b/rowSums(b)
    expFreq<- collapseSortMean(b2, numRep)
    save(expFreq, file=fileExpectedFreq)
}

plotLabels(plotName)
plotqq(expFreq, freqDataRef, paste0(plotName, " ", subName2) )


## dirichlet-multinomial mixture
for( m in which(whichIsP) ) {

    plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
    if(m==3){
        plotName<- gsub(" 1 components", "", plotName)
        plotName<- gsub("Mixture of ", "", plotName)
    }
    fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )
    if ( file.exists(fileExpectedFreq) && loadData ){
        load(fileExpectedFreq)
    } else{
        mModel<- maxModel[[m]]
        
        b <- rmdm(numRep*length(rowSumDataRef), rep(rowSumDataRef, numRep), mModel$f,
                        phi=mModel$params[,1], p=mModel$params[,c(2,4)])
        b2 <- b/rowSums(b)
        expFreq<- collapseSortMean(b2, numRep)
        save(expFreq, file=fileExpectedFreq)
    }
    plotName = gsub(" CHM1 Chr\\d\\d", "", plotName)
    plotName = gsub("Mixture of Dirichlet Multinomial (\\d) components", "Mixture of \\1 Dirichlet-Multinomials", plotName)
    plotName = gsub("Dirichlet Multinomial", "Dirichlet-Multinomial", plotName)

    main<- plotName
    plotLabels(plotName)
    plotqq(expFreq, freqDataRef, main)
}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")
} # match (p in 1:length(subNameList) )





# ####################################################################
# ##### Plot MS
# ####################################################################

# p<- 2

# subName<- subNameList[p]
# fullTitle<- fullTitleList[p]
# fullPath <- file.path(dataDir, subName)

# temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
# dataFull <- temp$dataFull
# #     dataRef <- temp$dataRef
# dataRefDirty <- temp$dataRefDirty

# maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)

# n <- rowSums(dataRefDirty)
# propRef <- dataRefDirty[,1]/n
# oo <- propRef > 0.8
# dataRefProp <- dataRefDirty[oo,]

# plotData<- dataRefDirty[,-2]
# plotData<- dataRefProp[,-2]
# rowSumDataRef<- rowSums(plotData)
# colSumDataRef<- colSums(plotData)

# freqDataRef<- plotData/rowSumDataRef
# subName2<- gsub("_C", " Chr", subName)


# whichIsP <- grepl("_[0-9]P",names(maxModel))
# header<- gsub("_" , "", names(maxModel) )
# header<- gsub("D" , " FD", header )
# header<- gsub("P" , " components", header ) #need this one
# header<-paste(subName2, header)


# plotTitle <- paste0(plotNamePrefix, "CHM1_MS.pdf") 
# qqplotFile<- file.path(figureDir, plotTitle)


# # pdf(file=qqplotFile, width=12, height=6, title=plotTitle)
# # par(mai=c(0.8,0.8,0.2,0.1), mfrow=c(1,2), cex.main=1.2^2,cex.lab=1.2, pty="s", omi=c(0,0,0,0) )

# ## simple multinomial (no bias)
# plotName<- "Multinomial"
# fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )

# if ( file.exists(fileExpectedFreq) && loadData ){
#     load(fileExpectedFreq)
# } else{
#     warning(paste("File doesn't exist:", fileExpectedFreq))
# }

# plotTitle <- paste0(plotNamePrefix, "CHM1_MS_",plotName,".png") 
# qqplotFile<- file.path(figureDir, plotTitle)

# png(file=qqplotFile, width = 600, height = 300, title=plotTitle)
# par(mai=c(0.8,0.8,0.2,0.1), mfrow=c(1,2), cex.main=1.2^2,cex.lab=1.2, pty="s", omi=c(0,0,0,0) )

# plotqq(expFreq, freqDataRef, "" )

# dev.off()

# ## dirichlet-multinomial mixture
# pName<-c("","","DM","","","2CMDM")
# for( m in which(whichIsP)[1:2] ) {

#     plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
#     if(m==3){
#         plotName<- gsub(" 1 components", "", plotName)
#         plotName<- gsub("Mixture of ", "", plotName)
#     }
#     fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )
#     if ( file.exists(fileExpectedFreq) && loadData ){
#         load(fileExpectedFreq)
#     } else{
#         warning(paste("File doesn't exist:", fileExpectedFreq))
#     }
    
# plotTitle <- paste0(plotNamePrefix, "CHM1_MS_",pName[m],".png") 
# qqplotFile<- file.path(figureDir, plotTitle)

# png(file=qqplotFile, width = 600, height = 300, title=plotTitle)
# par(mai=c(0.8,0.8,0.2,0.1), mfrow=c(1,2), cex.main=1.2^2,cex.lab=1.2, pty="s", omi=c(0,0,0,0) )

#     plotqq(expFreq, freqDataRef, "")
    
#     dev.off()
# }

# # dev.off()
# # embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")

