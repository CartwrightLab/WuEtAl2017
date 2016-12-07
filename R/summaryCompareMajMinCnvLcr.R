#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- TRUE
suppressPackageStartupMessages( source("./summaryFunctions.R") )
source("./summarySetup.R")


cnvFile1<- read.table(file.path(cnvDir, "Mills_deletions.csv"), sep=",", header=T)
cnvFile1<-cnvFile1[cnvFile1[,1]=="NA12878",]
cnvFile2<- read.table(file.path(cnvDir, "Mills_deletionsInsertions.csv"), sep=",", header=T)    
cnvFile2<-cnvFile2[cnvFile2[,1]=="NA12878",]


lc21<- read.table(paste0(lcDir, "lc_chr21.bed"), header=F)
lc10<- read.table(paste0(lcDir, "lc_chr10.bed"), header=F)    
lcRegionList<- list("21"=lc21, "10"=lc10)

################################################################################
## major vs minor
################################################################################
compareMajorMinorResultFileName<- paste0(latexTableDir, "Compare_MajMin_", projectName, ".tex")

#     Dataset & & FH & TH & FH percent & Not CNV & CNV & CNV percent & p-value & Not LCR & LCR & LCR percent & p-value \\\\ \\hline")

prefix<- paste0("\\begin{tabular}{|cc|cccc|cccc|}
    \\hline 
    Dataset & & Not CNV & CNV & CNV \\% & p-value & Not LCR & LCR & LCR \\% & p-value \\\\ \\hline")

sufix<- "\\end{tabular}"
    
cat(prefix, file=compareMajorMinorResultFileName, fill=TRUE)

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
    maxLikelihoodTable <- loadMaxLikelihoodTable(fullPath, subName, loadData, maxModel, dataFull, lowerLimit, upperLimit, isCEU, isRscriptMode)

    #CNV Mills et al.
    subset(cnvFile1, Chromosome == "chr10", select=c(SV_Start, 4))
    c1<- cnvFile1[cnvFile1[,2]==paste0("chr",chromosomeIndex), ]
    c2<- cnvFile2[cnvFile2[,2]==chromosomeIndex, ]

    cnvRegion<- c1[,3:4]
    cnvRegion<- rbind(cnvRegion, c2[,3:4])
    cnvRegionScale<- cnvRegion/1e6

    #LCR
    lcr<- lcRegionList[[chromosomeIndex]][,2:3]
    lcrScale<- lcr/1e6

    mlMaxIndex<- apply(maxLikelihoodTable[[BICIndexList[p]]], 1, which.max) 
    maxComp<- which.max(table(mlMaxIndex))
    predIndex<- (mlMaxIndex==maxComp)

    trueHetPos<- as.numeric(rownames(dataRef))
    potenHetPos<- as.numeric(rownames(dataRefDirty))
    predMajorPos<- potenHetPos[predIndex]
    predMinorPos<- potenHetPos[!predIndex]

    predMajorPosScale<- predMajorPos/1e6
    predMinorPosScale<- predMinorPos/1e6

    predMajorResult<- table(predMajorPos %in% trueHetPos)
    predMinorResult<-table(predMinorPos %in% trueHetPos)

    ##### CNV Major vs Minor
    cnvMajor<- sapply(predMajorPosScale, function(x){
        any(is.between(x, cnvRegionScale[,1], cnvRegionScale[,2]))
    })
    cnvMajorCountT<- sum(cnvMajor)
    cnvMajorCountF<- sum(!cnvMajor)

    cnvMinor<- sapply(predMinorPosScale, function(x){
        any(is.between(x, cnvRegionScale[,1], cnvRegionScale[,2]))
    })
    cnvMinorCountT<- sum(cnvMinor)
    cnvMinorCountF<- sum(!cnvMinor)

    combinedCnvCount<- matrix(c(table(cnvMajor), table(cnvMinor)), nrow=2, byrow=T)
    cnvResult<- fisher.test(combinedCnvCount)
    cnvStar<- ifelse(cnvResult$p.value< 1e-8, "*", " ") #1.46e-9 ... 0.0555



    ##### lcr Major vs Minor
    lcrMajor<- sapply(predMajorPosScale, function(x){
        any(is.between(x, lcrScale[,1], lcrScale[,2]))
    })
    lcrMajorCountT<- sum(lcrMajor)
    lcrMajorCountF<- sum(!lcrMajor)

    lcrMinor<- sapply(predMinorPosScale, function(x){
        any(is.between(x, lcrScale[,1], lcrScale[,2]))
    })
    lcrMinorCountT<- sum(lcrMinor)
    lcrMinorCountF<- sum(!lcrMinor)

    combinedLcrCount<- matrix(c(table(lcrMajor), table(lcrMinor)), nrow=2, byrow=T)
    lcrResult<- fisher.test(combinedLcrCount)
    lcrStar<- ifelse(lcrResult$p.value< 1e-6, "*", " ")  #1.64e-7 1.93e-33



    cnvStar<- ""
    lcrStar<- ""    
    formatString<- paste0(
        paste(fullTitle, "Major", #paste(predMajorResult, collapse=" & "), 
#             formatC(100*predMajorResult[1]/sum(predMajorResult), digits=1 ,format="f", flag="#"), 
            cnvMajorCountF, cnvMajorCountT, 
            paste0(formatC(100*cnvMajorCountT/length(cnvMajor), digits=1, format="f", flag="#"), cnvStar),
    #         formatC(cnvMajorCountT/length(cnvMajor), digits=4, format="f", flag="#"),
            format(cnvResult$p.value, digits=3, format="f", flag="#"),
            lcrMajorCountF, lcrMajorCountT, 
            paste0(formatC(100*lcrMajorCountT/length(lcrMajor), digits=1, format="f", flag="#"), lcrStar),
    #         formatC(lcrMajorCountT/length(lcrMajor), digits=3, format="f", flag="#"),
            format(lcrResult$p.value, digits=3, format="f", flag="#"),
            sep=" & "),
        "\\\\ ")
    cat(formatString, file=compareMajorMinorResultFileName, append=TRUE, fill=TRUE)
    
    formatString<- paste0(
        paste("", "Minor", #paste(predMinorResult, collapse=" & "), 
#             formatC(100*predMinorResult[1]/sum(predMinorResult), digits=1 ,format="f", flag="#"), 
            cnvMinorCountF, cnvMinorCountT, 
            formatC(100*cnvMinorCountT/length(cnvMinor), digits=1, format="f", flag="#"), 
            "",
            lcrMinorCountF, lcrMinorCountT, 
            formatC(100*lcrMinorCountT/length(lcrMinor), digits=1, format="f", flag="#"),
            "",
            sep=" & "),
        " \\\\ \\hline")
    cat(formatString, file=compareMajorMinorResultFileName, append=TRUE, fill=TRUE)

}

cat(sufix, file=compareMajorMinorResultFileName, fill=T, append=T)
