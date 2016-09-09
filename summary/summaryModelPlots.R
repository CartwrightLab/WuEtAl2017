#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- TRUE
suppressPackageStartupMessages( source("./summaryFunctions.R") )
source("./summarySetup.R")

msTable<- vector(mode="list", length=length(subNameList)*2)
AicThIndex <- c(4,3,3,4,4,4,3,4)
BicThIndex <- BICIndexList/2

AicPhIndex <- c(5,6,6,6,6,6,3,6)
BicPhIndex <- c(3,4,4,5,4,5,3,6)

allModel2<- vector(length=length(subNameList), mode="list")
for(p in length(subNameList):1 ){

    subName<- subNameList[p]
    fullTitle<- fullTitleList[p]
    fullPath <- file.path(dataDir, subName)

    temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
    dataFull <- temp$dataFull
    dataRef <- temp$dataRef
    dataRefDirty <- temp$dataRefDirty

    
    maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)
    

    allModel2[[p]]<- 
    sapply(maxModel[3], function(x){
        sd<- sqrt(c( x$covar["phi.1", "phi.1"], x$covar["phi.2", "phi.2"] ))
        r<- cbind("phi"=x$params[,1:1], sd, prop=x[["f"]] )
        newOrder<-rev(order(r[,3]))
        r<- r[newOrder,]
        return(r)
    }, simplify=F )
}
    
allResult<- lapply(allModel2, function(x){x[[1]]})
names(allResult)<- subNameList

plotData<- sapply(allResult, function(x){
                    m<- x[1,1]
                    sd<- x[1,2]
                    r<- c(mean=m, interval=1.96*sd, lower=m-1.96*sd, upper=m+1.96*sd)
        })

newOrder<- c( rev(grep("C21", colnames(plotData))), rev(grep("C10", colnames(plotData)))  )
plotData2<- plotData[, newOrder]

require(gplots)
phiPlotFile<- file.path(figureDir, paste0("phiPlot.pdf") )
pdf(file=phiPlotFile, width=6, height=6, title=phiPlotFile,point) 
par(mar=c(7,3,0,0), mgp=c(2,1,0), cex=1)

plotCI(x=plotData2["mean",], uiw=plotData2["interval",], minbar=0, ylab=expression(phi), xaxt="n", xlab="") #xlab="Dataset")
mtext("Dataset", 1, line=6)
axis(1,1:8, lab=gsub("_", " ", colnames(plotData2)), las=2)

dev.off()


