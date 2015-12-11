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

## lc(r) - low complexity (region)

# source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/mdm.R")

isCEU <- FALSE
isCEU <- TRUE

dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10
loadData<- TRUE


latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

lcDir <- "/home/steven/Postdoc2/Project_MDM/LowComp/"

if(isCEU){
    pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
    subNameList<- c(
    "CEU10_C10", "CEU10_C21", 
    "CEU11_C10", "CEU11_C21",
    "CEU12_C10", "CEU12_C21",
    "CEU13_C10", "CEU13_C21"
    )

    fullTitleList<- c(
    "CEU10 Chr10", "CEU10 Chr21",
    "CEU11 Chr10", "CEU11 Chr21",
    "CEU12 Chr10", "CEU12 Chr21",
    "CEU13 Chr10", "CEU13 Chr21"
    )
    projectName<- "CEU"
} else {
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
    projectName<- "CHM1"
}
setwd(pwd)


lc21<- read.table(paste0(lcDir, "lc_chr21.bed"), header=F)
lc10<- read.table(paste0(lcDir, "lc_chr10.bed"), header=F)    
lcRegionList<- list("21"=lc21, "10"=lc10)

is.between <- function(x, a, b) {
    ( (x - a)  *  (b - x) ) >= 0
}

################################################################################
##### TH vs FH
################################################################################

lcrResultFileName<- paste0(latexDir, "LowComp_", projectName, ".tex")
prefix<- paste0("\\begin{tabular}{|cc|ccc|c|}
    \\hline 
    Dataset & & Not CNV & CNV & CNV proportion & p-value \\\\ \\hline")
sufix<- "\\end{tabular}"
    

cat(prefix, file=lcResultFileName, fill=TRUE)
for(p in length(subNameList):1 ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
cat(subName, fill=T)

lcr<- lcRegionList[[chromosomeIndex]][,2:3]
lcrScale<- lcr/1e6

if(isCEU){
    hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
    dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
    dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
    dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE)
} else{
    hets_byref<- file.path("base_count_meta_subsample") 
    dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
    # dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU=isCEU)
    dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU=isCEU)
}

## True H vs False H
trueHetName<- as.numeric(rownames(dataRef))
potHetName<- as.numeric(rownames(dataRefDirty))
falsePosIndex <- which(! potHetName %in% trueHetName)
falsePosPosition<- as.integer(potHetName[falsePosIndex])

falsePosPositionScale<- falsePosPosition/1e6
trueHetPositionScale<- trueHetName/1e6

# falseTF<- sapply(falsePosPosition, function(x){
#         any((x >= lcr[,1] & x <= lcr[,2]  ))  

falseTF<- sapply(falsePosPositionScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
    
falseTFCountT<- sum(falseTF)
falseTFCountF<- sum(!falseTF)
# summary(falseTF)
# cat(falseTFCountF,  falseTFCountT, formatC(falseTFCountT/length(falseTF)), "\n")

# trueTF<- sapply(trueHetName, function(x){
#         if( any((x >= lcr[,1] & x <= lcr[,2]  ))  ){
#             return(TRUE)
#         }
#         return(FALSE)
#     })
# })

trueTF<- sapply(trueHetPositionScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
trueTFcountT<- sum(trueTF)
trueTFcountF<- sum(!trueTF)
# summary(trueTF)
# cat(trueTFcountF,  trueTFcountT, formatC(trueTFcountT/length(trueTF)), "\n")
 
combinedLcrCount<- matrix(c(table(falseTF), table(trueTF)), nrow=2, byrow=T)
r<- fisher.test(combinedLcrCount)

formatString<- paste0(
    paste(fullTitle, "FH", falseTFCountF,  falseTFCountT, 
        formatC(falseTFCountT/length(falseTF), digits=3), 
        formatC(r$p.value, digits=3), sep=" & "), 
    " \\\\")
cat(formatString, file=lcrResultFileName, append=TRUE, fill=TRUE)
    
formatString<- paste0(
    paste("", "TH", trueTFcountF,  trueTFcountT, 
        formatC(trueTFcountT/length(trueTF), digits=3), "", sep=" & "),
    " \\\\ \\hline") 
cat(formatString, file=lcrResultFileName, append=TRUE, fill=TRUE)

}

cat(sufix, file=lcrResultFileName, fill=T, append=T)


################################################################################
## major vs minor
################################################################################
lcMajorMinorResultFileName<- paste0(latexDir, "LowComp_MajorMinorCat_", projectName, ".tex")

prefix<- paste0("\\begin{tabular}{|c|ccc|ccc|c|}
    \\hline 
    Dataset & FH & TH & FH proportion & Not CNV & CNV & CNV proportion & p-value \\\\ \\hline")

sufix<- "\\end{tabular}"
    

cat(prefix, file=lcMajorMinorResultFileName, fill=TRUE)
for(p in length(subNameList):1 ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
cat(subName, fill=T)

lcr<- lcRegionList[[chromosomeIndex]][,2:3]
lcrScale<- lcr/1e6

if(isCEU){
    hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
    dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
    dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
    dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE)
} else{
    hets_byref<- file.path("base_count_meta_subsample") 
    dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
    # dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU=isCEU)
    dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU=isCEU)
}

# setwd(fullPath)
maxModel<- extractMaxModel(fullPath)
fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTableFull.RData")
if ( file.exists(fileMaxLikelihoodTabel) && loadData ){
    load(fileMaxLikelihoodTabel)
} 

mlMaxIndex<- apply(maxLikelihoodTable[[BICIndex]], 1, which.max) 
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

##### Major
majorTF<- sapply(predMajorPosScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
    
# summary(majorTF)
majorCount_T<- sum(majorTF)
majorCount_F<- sum(!majorTF)
# cat(table(majorTF), formatC(majorCount_T/length(majorTF)), "\n")

#### Minor

minorTF<- sapply(predMinorPosScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
# table(minorTF)
MinorCount_T<- sum(minorTF)
MinorCount_F<- sum(!minorTF)
# cat(table(minorTF), formatC(MinorCount_T/length(minorTF)), "\n")

combinedCnvCount<- matrix(c(table(majorTF), table(minorTF)), nrow=2, byrow=T)
r<- fisher.test(combinedCnvCount)

formatString<- paste0(
    paste(fullTitle, paste(predMajorResult, collapse=" & "), 
        formatC(predMajorResult[1]/sum(predMajorResult), digits=3), 
        majorCount_F,  majorCount_T, formatC(majorCount_T/length(majorTF), digits=3),
        formatC(r$p.value, digits=3), sep=" & "),
    "\\\\ ")
cat(formatString, file=lcMajorMinorResultFileName, append=TRUE, fill=TRUE)
    
formatString<- paste0(
    paste(" ", paste(predMinorResult, collapse=" & "), 
        formatC(predMinorResult[1]/sum(predMinorResult), digits=3), 
        MinorCount_F,  MinorCount_T, formatC(MinorCount_T/length(minorTF), digits=3), " ", sep=" & "),
    " \\\\ \\hline")
cat(formatString, file=lcMajorMinorResultFileName, append=TRUE, fill=TRUE)

}

cat(sufix, file=lcMajorMinorResultFileName, fill=T, append=T)
