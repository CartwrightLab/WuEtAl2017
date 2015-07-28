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
isCEU <- TRUE

dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10
loadData<- TRUE


latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

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
    "CEU 2010 Chromosome 10", "CEU 2010 Chromosome 21",
    "CEU 2011 Chromosome 10", "CEU 2011 Chromosome 21",
    "CEU 2012 Chromosome 10", "CEU 2012 Chromosome 21",
    "CEU 2013 Chromosome 10", "CEU 2013 Chromosome 21"
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


p<- 8
p<- 1

cnvResultFileName<- paste0(latexDir, "CNV_", projectName, ".tex")

prefix<- paste0("\\begin{tabular}{|c|c|c|c|}
    \\hline \\multicolumn{4}{|c|}{Copy number variation (CNV) summary}\\\\ \\hline
    Dataset & Not CNV & CNV & CNV Percentage \\\\ \\hline")
sufix<- "\\end{tabular}"

    
cat(prefix, file=cnvResultFileName, fill=TRUE)
for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
cnvFile<- read.table(paste0(cnvDir, "CNV_C", chromosomeIndex, "_list"))

# cnvRegion<- cnvFile[, 2:3]
# table(cnvFile[,1])
cnvRegion<- cnvFile[cnvFile[,1]=="copy_number_variation" | 
                    cnvFile[,1]=="copy_number_gain" | 
                    cnvFile[,1]=="copy_number_loss", 2:3]

setwd(fullPath)
maxModel<- extractMaxModel(fullPath)

hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE)

# cnvFile[,1]



trueHetName<- as.numeric(rownames(dataRef))
potHetName<- as.numeric(rownames(dataRefDirty))
falsePosIndex <- which(! potHetName %in% trueHetName)
falsePosPosition<- as.integer(potHetName[falsePosIndex])
uniqueCnvRegion<- unique(cnvRegion)

cnvTF<- vector(length=length(falsePosPosition))
for(i in 1:length(falsePosPosition) ){
    index<- falsePosPosition[i]
    if( any(which(index >= uniqueCnvRegion[,1] & index <= uniqueCnvRegion[,2]  ))  ){
        cnvTF[i] <- TRUE
    }
}

summary(cnvTF)
count_T<- sum(cnvTF)
count_F<- length(cnvTF)-count_T
cat(paste0(fullTitle, " False positive sites"), " & ",
    paste(count_F,  count_T, formatC(count_T/length(cnvTF)), sep=" & "),
    " \\\\ " , file=cnvResultFileName, append=TRUE, fill=TRUE)


trueHetName<- as.numeric(trueHetName)
allTF<- vector(length=length(trueHetName))
for(i in 1:length(trueHetName) ){
    index<- trueHetName[i]
    if( any(which(index >= uniqueCnvRegion[,1] & index <= uniqueCnvRegion[,2]  ))  ){
        cnvTF[i] <- TRUE
    }
}
summary(allTF)
count_T<- sum(allTF)
count_F<- length(allTF)-count_T

cat(paste0(fullTitle, " All heterozygous"), " & ",
    paste(count_F,  count_T, formatC(count_T/length(cnvTF)), sep=" & "),
    " \\\\ \\hline" , file=cnvResultFileName, append=TRUE, fill=TRUE)
# unique( cbind( rep(1:2, 60), rep(1:5, 24)) )

}
cat(sufix, file=cnvResultFileName, fill=T, append=T)
