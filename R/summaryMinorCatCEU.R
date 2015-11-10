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

cnvResultFileName<- paste0(latexDir, "MapMinorCat_", projectName, ".tex")
BICIndexList<- c(2,2,2,2,2,2,2,3)*2 #TH
# cnvResultFileName<- paste0(latexDir, "MapMinorCatBICPH_", projectName, ".tex")
# BICIndexList<- c(3,4,4,5,4,5,3,6)*2 #PH

prefix<- paste0("\\begin{tabular}{|c|ccc|ccc|}
    \\hline \\multicolumn{7}{|c|}{Non-major categories summary}\\\\ \\hline
    Dataset & FH & TH & FH proportion & CNV & CNV & CNV proportion \\\\ \\hline")
sufix<- "\\end{tabular}"


# cnvFile<- read.table(paste0(cnvDir, "DGV_GRCh37_hg19_variants_subset.txt"), header=T)
        
cnvFile1<- read.table(paste0(cnvDir, "Mills_deletions.csv"), sep=",", header=T)
cnvFile1<-cnvFile1[cnvFile1[,1]=="NA12878",]
cnvFile2<- read.table(paste0(cnvDir, "Mills_deletionsInsertions.csv"), sep=",", header=T)    
cnvFile2<-cnvFile2[cnvFile2[,1]=="NA12878",]

cat(prefix, file=cnvResultFileName, fill=TRUE)


for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
BICIndex<- BICIndexList[p]

# cnvFile<- read.table(paste0(cnvDir, "CNV_C", chromosomeIndex, "_list"))


#Mills et al.
subset(cnvFile1, Chromosome == "chr10", select=c(SV_Start, 4))
c1<- cnvFile1[cnvFile1[,2]==paste0("chr",chromosomeIndex), ]
c2<- cnvFile2[cnvFile2[,2]==chromosomeIndex, ]

cnvRegion<- c1[,3:4]
cnvRegion<- rbind(cnvRegion, c2[,3:4])
uniqueCnvRegion<- cnvRegion


# cnvRegion<- cnvFile[, 2:3]
# table(cnvFile[,1])
# cnvRegion<- cnvFile[cnvFile[,1]=="copy_number_variation" | 
#                     cnvFile[,1]=="copy_number_gain" | 
#                     cnvFile[,1]=="copy_number_loss", 2:3]

#DGV
# cnvChromosome<- cnvFile[cnvFile[,2]==chromosomeIndex & cnvFile[,5]=="CNV",]
# cnvRegion<- cnvChromosome[,3:4]
# uniqueCnvRegion<- unique(cnvRegion)

setwd(fullPath)
maxModel<- extractMaxModel(fullPath, isCEU=isCEU)


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

# maxLikelihoodList<- maxLikelihoodTable
# result<- sapply(maxLikelihoodList, function(x){
# #                 prop <- table(apply(x,1,which.max))
#                 prop <- tabulate( apply(x,1,which.max) , nbins=NCOL(x) )
#                 return(prop/sum(prop))
#         })
# 

fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTableFull.RData")
if ( file.exists(fileMaxLikelihoodTabel) && loadData ){
    load(fileMaxLikelihoodTabel)
} 
# else{
#     maxLikelihoodTable<- calculateEachLikelihood(maxModel, dataFull, lowerLimit=lowerLimit, upperLimit=upperLimit, isCEU=isCEU)#, numData=100)
#     attr(maxLikelihoodTable, "title")<- subName
#     save(maxLikelihoodTable, file=file.path(fullPath, "maxLikelihoodTableFull.RData"))
# }

# mlMaxIndex<- sapply(maxLikelihoodTable, function(x){
#                 prop <-  apply(x,1,which.max) 
#             })
mlMaxIndex<- apply(maxLikelihoodTable[[BICIndex]], 1, which.max) 


maxComp<- which.max(table(mlMaxIndex))
predIndex<- (mlMaxIndex==maxComp)

table(predIndex)

trueHetPos<- as.numeric(rownames(dataRef))
potenHetPos<- as.numeric(rownames(dataRefDirty))
predMajorPos<- potenHetPos[predIndex]
predMinorPos<- potenHetPos[!predIndex]

predResult<-table(predMinorPos %in% trueHetPos)

cnvTF<- sapply(predMinorPos, function(x){
            if( any((x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  ))  ){
                return(TRUE)
            }
            return(FALSE)
        })


table(cnvTF)
count_T<- sum(cnvTF)
count_F<- sum(!cnvTF)
cat(table(cnvTF), formatC(count_T/length(cnvTF)), "\n")


cat(paste0(fullTitle, " & "),
    paste(predResult, " & ", collapse=""), formatC(predResult[1]/sum(predResult)), " & ",
    paste(count_F,  count_T, formatC(count_T/length(cnvTF)), sep=" & "),
    " \\\\ \\hline" , file=cnvResultFileName, append=TRUE, fill=TRUE)

    
##########
##########

allTF<- sapply(predMajorPos, function(x){
        if( any((x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  ))  ){
            return(TRUE)
        }
        return(FALSE)
    })
    
summary(allTF)
count_T<- sum(allTF)
count_F<- sum(!allTF)
cat(table(allTF), formatC(count_T/length(allTF)), "\n")

predMajorResult<- table(predMajorPos %in% trueHetPos)

cat(paste0(fullTitle, " Major component"), " & ",
    paste(predMajorResult, " & ", collapse=""), 
    formatC(predMajorResult[1]/sum(predMajorResult)), " & ",
    paste(count_F,  count_T, formatC(count_T/length(allTF)), sep=" & "),
    " \\\\ \\hline" , file=cnvResultFileName, append=TRUE, fill=TRUE)

    
# cat(paste0(fullTitle, " Major component"), " & ",
#     paste(count_F,  count_T, formatC(count_T/length(cnvTF)), sep=" & "),
#     " \\\\ \\hline" , file=cnvResultFileName, append=TRUE, fill=TRUE)


}
cat(sufix, file=cnvResultFileName, fill=T, append=T)


