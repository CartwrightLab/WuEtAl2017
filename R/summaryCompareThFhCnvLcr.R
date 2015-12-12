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
##### TH vs FH
################################################################################

ThFhResultFileName<- paste0(latexDir, "Compare_ThFh_", projectName, ".tex")
prefix<- paste0("\\begin{tabular}{|cc|cccc|cccc|}
    \\hline
    Dataset & & Not CNV & CNV & CNV proportion & p-value & Not LCR & LCR & LCR proportion & p-value\\\\ \\hline")
sufix<- "\\end{tabular}"
    

cat(prefix, file=ThFhResultFileName, fill=TRUE)
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

## True H vs False H
trueHetName<- as.numeric(rownames(dataRef))
potHetName<- as.numeric(rownames(dataRefDirty))
falsePosIndex <- which(! potHetName %in% trueHetName)
falsePosPosition<- as.integer(potHetName[falsePosIndex])

falsePosPositionScale<- falsePosPosition/1e6
trueHetPositionScale<- trueHetName/1e6


#CNV
cnvFalseH<- sapply(falsePosPositionScale, function(x){
    any(is.between(x, cnvRegionScale[,1], cnvRegionScale[,2] ))
})
cnvFalseHCountT<- sum(cnvFalseH)
cnvFalseHCountF<- sum(!cnvFalseH)

cnvTrueH<- sapply(trueHetPositionScale, function(x){
    any(is.between(x, cnvRegionScale[,1], cnvRegionScale[,2] ))
})
cnvTrueHCountT<- sum(cnvTrueH)
cnvTrueHCountF<- sum(!cnvTrueH)

#LCR
lcrFalseH<- sapply(falsePosPositionScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
lcrFalseHCountT<- sum(lcrFalseH)
lcrFalseHCountF<- sum(!lcrFalseH)


lcrTrueH<- sapply(trueHetPositionScale, function(x){
    any(is.between(x, lcrScale[,1], lcrScale[,2]))
})
lcrTrueHCountT<- sum(lcrTrueH)
lcrTrueHCountF<- sum(!lcrTrueH)
 

combinedCnvCount<- matrix(c(table(cnvFalseH), table(cnvTrueH)), nrow=2, byrow=T)
cnvResult<- fisher.test(combinedCnvCount)
cnvStar<- ifelse(cnvResult$p.value< 1e-8, "*", " ") #1.46e-9 ... 0.0555

combinedLcrCount<- matrix(c(table(lcrFalseH), table(lcrTrueH)), nrow=2, byrow=T)
lcrResult<- fisher.test(combinedLcrCount)
lcrStar<- ifelse(lcrResult$p.value< 1e-6, "*", " ")  #1.64e-7 1.93e-33

cnvStar<- ""
lcrStar<- ""
formatString<- paste0(
    paste(fullTitle, "TH", 
        cnvTrueHCountF,  cnvTrueHCountT, 
        paste0(formatC(cnvTrueHCountT/length(cnvTrueH), digits=4, format="f", flag="#"), cnvStar),
        format(cnvResult$p.value, digits=3, format="f", flag="#"),
        lcrTrueHCountF,  lcrTrueHCountT, 
        paste0(formatC(lcrTrueHCountT/length(lcrTrueH), digits=3, format="f", flag="#"), lcrStar),
        format(lcrResult$p.value, digits=3, format="f", flag="#"),
        sep=" & "),
    " \\\\")
cat(formatString, file=ThFhResultFileName, append=TRUE, fill=TRUE)
    
formatString<- paste0(
    paste("", "FH",
        cnvFalseHCountF,  cnvFalseHCountT, 
        formatC(cnvFalseHCountT/length(cnvFalseH), digits=4, format="f", flag="#"),
        "",
        lcrFalseHCountF,  lcrFalseHCountT, 
        formatC(lcrFalseHCountT/length(lcrFalseH), digits=3, format="f", flag="#"),
        "",
        sep=" & "),
    " \\\\ \\hline") 
cat(formatString, file=ThFhResultFileName, append=TRUE, fill=TRUE)

}

cat(sufix, file=ThFhResultFileName, fill=T, append=T)
