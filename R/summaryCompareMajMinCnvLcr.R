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
cnvDir <- "/home/steven/Postdoc2/Project_MDM/CNV/"

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


BICIndexList<- c(2,2,2,2,2,2,2,3)*2 #TH

cnvFile1<- read.table(paste0(cnvDir, "Mills_deletions.csv"), sep=",", header=T)
cnvFile1<-cnvFile1[cnvFile1[,1]=="NA12878",]
cnvFile2<- read.table(paste0(cnvDir, "Mills_deletionsInsertions.csv"), sep=",", header=T)    
cnvFile2<-cnvFile2[cnvFile2[,1]=="NA12878",]


lc21<- read.table(paste0(lcDir, "lc_chr21.bed"), header=F)
lc10<- read.table(paste0(lcDir, "lc_chr10.bed"), header=F)    
lcRegionList<- list("21"=lc21, "10"=lc10)

is.between <- function(x, a, b) {
    ( (x - a)  *  (b - x) ) >= 0
}

################################################################################
## major vs minor
################################################################################
compareMajorMinorResultFileName<- paste0(latexDir, "Compare_MajMin_", projectName, ".tex")

prefix<- paste0("\\begin{tabular}{|cc|ccc|cccc|cccc|}
    \\hline 
    Dataset & & FH & TH & FH proportion & Not CNV & CNV & CNV proportion & p-value & Not LCR & LCR & LCR proportion & p-value \\\\ \\hline")

sufix<- "\\end{tabular}"
    
cat(prefix, file=compareMajorMinorResultFileName, fill=TRUE)

for(p in length(subNameList):1 ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
cat(subName, fill=T)

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
    paste(fullTitle, "Major", paste(predMajorResult, collapse=" & "), 
        formatC(predMajorResult[1]/sum(predMajorResult), digits=3 ,format="f", flag="#"), 
        cnvMajorCountF, cnvMajorCountT, 
        paste0(formatC(cnvMajorCountT/length(cnvMajor), digits=4, format="f", flag="#"), cnvStar),
#         formatC(cnvMajorCountT/length(cnvMajor), digits=4, format="f", flag="#"),
        format(cnvResult$p.value, digits=3, format="f", flag="#"),
        lcrMajorCountF, lcrMajorCountT, 
        paste0(formatC(lcrMajorCountT/length(lcrMajor), digits=3, format="f", flag="#"), lcrStar),
#         formatC(lcrMajorCountT/length(lcrMajor), digits=3, format="f", flag="#"),
        format(lcrResult$p.value, digits=3, format="f", flag="#"),
        sep=" & "),
    "\\\\ ")
cat(formatString, file=compareMajorMinorResultFileName, append=TRUE, fill=TRUE)
    
formatString<- paste0(
    paste("", "Minor", paste(predMinorResult, collapse=" & "), 
        formatC(predMinorResult[1]/sum(predMinorResult), digits=3 ,format="f", flag="#"), 
        cnvMinorCountF, cnvMinorCountT, 
        formatC(cnvMinorCountT/length(cnvMinor), digits=3, format="f", flag="#"), 
        "",
        lcrMinorCountF, lcrMinorCountT, 
        formatC(lcrMinorCountT/length(lcrMinor), digits=3, format="f", flag="#"),
        "",
        sep=" & "),
    " \\\\ \\hline")
cat(formatString, file=compareMajorMinorResultFileName, append=TRUE, fill=TRUE)

}

cat(sufix, file=compareMajorMinorResultFileName, fill=T, append=T)
