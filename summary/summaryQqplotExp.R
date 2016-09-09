#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- TRUE
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
   
    chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)

    temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
    dataFull <- temp$dataFull
    dataRef <- temp$dataRef
    dataRefDirty <- temp$dataRefDirty

    maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)

    ## QQ plots
    rowSumDataRef<- rowSums(dataRef)
    colSumDataRef<- colSums(dataRef)
    rowSumDataRefDirty<- rowSums(dataRefDirty)

    freqDataRef<- dataRef/rowSumDataRef
    freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty

    subName2<- gsub("_C", " Chr", subName)


    whichIsDirty <- grepl("_[0-9]D",names(maxModel))
    header<- gsub("hets_" , "", names(maxModel) )
    header<- gsub(".*_", "", header)
    header<- gsub("$", " components", header)
    header<-paste(subName2, header)

    qqplotFile<- file.path(figureDir, paste0(plotNamePrefix, subName, ".pdf") )

    mtop = 0.4
    mbottom = 0.3
    mleft = 0.4
    mright = 0
    width = 7.5
    sz = (width-mleft-mright)/3
    height = mtop+mbottom+sz

    m = c(
        0,1,1,1,0,
        2,4,5,6,0,
        0,3,3,3,0
    );
    m = matrix(m,3,5,byrow=TRUE)
    w = c(mleft,sz,sz,sz,mright)
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
        p1<- (colSumDataRef[1]+colSumDataRef[2])/2
        prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
        ll <- sum(log(prob)*colSumDataRef)
        # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
        b <- rmultinomial(length(rowSumDataRef)*numRep, rep(rowSumDataRef, numRep), prob)
        b2 <- b/rowSums(b)
        expFreq<- collapseSortMean(b2, numRep)
        save(expFreq, file=fileExpectedFreq)
    }

    plotLabels(plotName)
    plotqq(expFreq, freqDataRef, paste0(plotName, " ", subName2) )


    ## multinomial (with ref bias)
    plotName<- "Biased Multinomial"
    fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )
    
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

    plotLabels(plotName)
    plotqq(expFreq, freqDataRef, paste0(plotName, " ", subName2) )


    ## dirichlet-multinomial mixture
    for( m in which(! whichIsDirty) ) {
        
        plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
        if(m==1){
            plotName<- gsub(" 1 components", "", plotName)
            plotName<- gsub("Mixture of ", "", plotName)
        }

        fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )

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
        plotName = gsub(" CEU\\d\\d Chr\\d\\d", "", plotName)
        plotName = gsub("Mixture of Dirichlet Multinomial (\\d) components", "Mixture of \\1 Dirichlet-Multinomials", plotName)
        plotName = gsub("Dirichlet Multinomial", "Dirichlet-Multinomial", plotName)

        plotLabels(plotName)
        plotqq(expFreq, freqDataRef, plotName)
    }


    dev.off()
    embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")
    
} ## match (p in 1:length(subNameList) )


# ####################################################################
# ##### Plot MS
# ####################################################################

# # qqplotFile<- file.path(figureDir, paste0(plotNamePrefix, "CEU_MS.pdf") )
# # pdf(file=qqplotFile, width=15, height=6, title=qqplotFile)

# # qqplotFile<- file.path(figureDir, paste0(plotNamePrefix, "CEU_MS.png") )
# # png(filename = qqplotFile, width = 900, height = 300)
# # par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), cex.main=1.2^4,cex.lab=1.2^2, pty="s", omi=c(0,0,0,0) )
    
# p<- 8

#     subName<- subNameList[p]
# #     fullTitle<- fullTitleList[p]
#     fullPath <- file.path(dataDir, subName)
   
# #     chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)

#     temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
# #     dataFull <- temp$dataFull
#     dataRef <- temp$dataRef
# #     dataRefDirty <- temp$dataRefDirty

#     maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)

# #     ## QQ plots
#     rowSumDataRef<- rowSums(dataRef)
# #     colSumDataRef<- colSums(dataRef)
# #     rowSumDataRefDirty<- rowSums(dataRefDirty)

#     freqDataRef<- dataRef/rowSumDataRef
# #     freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty

#     subName2<- gsub("_C", " Chr", subName)


#     whichIsDirty <- grepl("_[0-9]D",names(maxModel))
#     header<- gsub("hets_" , "", names(maxModel) )
#     header<- gsub(".*_", "", header)
#     header<- gsub("$", " components", header)
#     header<-paste(subName2, header)

# ## simple multinomial (no bias)
# plotName<- "Multinomial"
# fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )

# if ( file.exists(fileExpectedFreq) && loadData ){
#     load(fileExpectedFreq)
# } else{
#     warning(paste("File doesn't exist:", fileExpectedFreq))
# }
# qqplotFile<- file.path(figureDir, paste0(plotNamePrefix, "CEU_MS_",plotName,".png") )
# png(filename = qqplotFile, width = 900, height = 300)
# par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), cex.main=1.2^4,cex.lab=1.2^2, pty="s", omi=c(0,0,0,0) )

# plotqq(expFreq, freqDataRef, "")

# dev.off()

# ## dirichlet-multinomial mixture
# pName<- c("DM","","","","3CMDM")
# for( m in which(! whichIsDirty)[c(1,3)] ) {
    
#     plotName<- paste0("Mixture of Dirichlet Multinomial ", header[m])
#     if(m==1){
#         plotName<- gsub(" 1 components", "", plotName)
#         plotName<- gsub("Mixture of ", "", plotName)
#     }

    
#     fileExpectedFreq <- file.path(fullPath, paste0(subName, "_", gsub(" ", "_", plotName), "_expectedFreq.RData") )
#     if ( file.exists(fileExpectedFreq) && loadData ){
#         load(fileExpectedFreq)
#     } else{
#         warning(paste("File doesn't exist:", fileExpectedFreq))
#     }
#     qqplotFile<- file.path(figureDir, paste0(plotNamePrefix, "CEU_MS_",pName[m],".png") )
#     png(filename = qqplotFile, width = 900, height = 300)
#     par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), cex.main=1.2^4,cex.lab=1.2^2, pty="s", omi=c(0,0,0,0) )

#     plotqq(expFreq, freqDataRef, "")
    
#     dev.off()
# }


# dev.off()
# # embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")

