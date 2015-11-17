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


p<- 8
p<- 1

cnvResultFileName<- paste0(latexDir, "CNV_", projectName, ".tex")

prefix<- paste0("\\begin{tabular}{|cc|ccc|c|}
    \\hline 
    Dataset & & Not CNV & CNV & CNV proportion & p-value \\\\ \\hline")
sufix<- "\\end{tabular}"

# cnvFile<- read.table(paste0(cnvDir, "DGV_GRCh37_hg19_variants_subset.txt"), header=T) #all human
# cnvFile<- read.table(paste0(cnvDir, "NA12878.wgs.illumina_platinum.20140404.svs_v2.vcf"))
    
cnvFile1<- read.table(paste0(cnvDir, "Mills_deletions.csv"), sep=",", header=T)
cnvFile1<-cnvFile1[cnvFile1[,1]=="NA12878",]
cnvFile2<- read.table(paste0(cnvDir, "Mills_deletionsInsertions.csv"), sep=",", header=T)    
cnvFile2<-cnvFile2[cnvFile2[,1]=="NA12878",]

cat(prefix, file=cnvResultFileName, fill=TRUE)
for(p in length(subNameList):1 ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)

# cnvFile<- read.table(paste0(cnvDir, "CNV_C", chromosomeIndex, "_list"))
# cnvRegion<- cnvFile[, 2:3]
# table(cnvFile[,1])
#  cnvRegion<- cnvFile[cnvFile[,1]=="copy_number_variation" |
#                      cnvFile[,1]=="copy_number_gain" | 
#                     cnvFile[,1]=="copy_number_loss"
#                      , 2:3]

#DGV_GRCh37_hg19_variants_subset.txt
# cnvChromosome<- cnvFile[cnvFile[,2]==chromosomeIndex & cnvFile[,5]=="CNV",]
# cnvRegion<- cnvChromosome[,3:4]
# uniqueCnvRegion<- unique(cnvRegion)


# cnvFile<- read.table(paste0(cnvDir, "NA12878_C", chromosomeIndex))
# text<-as.matrix(cnvFile[,8])
# endPos<- apply(text,1, function(x){
#     x<-as.vector(x)
#     as.numeric( gsub("caller.*END=", "", x) )
# })
# 
# cnvRegion<- cbind(cnvFile[, 2], endPos)
# uniqueCnvRegion<-cnvRegion



#Mills et al.
subset(cnvFile1, Chromosome == "chr10", select=c(SV_Start, 4))
c1<- cnvFile1[cnvFile1[,2]==paste0("chr",chromosomeIndex), ]
c2<- cnvFile2[cnvFile2[,2]==chromosomeIndex, ]

cnvRegion<- c1[,3:4]
cnvRegion<- rbind(cnvRegion, c2[,3:4])
uniqueCnvRegion<- cnvRegion






# dd<- apply(cnvRegion,1,diff)

setwd(fullPath)
maxModel<- extractMaxModel(fullPath)


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

# cnvFile[,1]



trueHetName<- as.numeric(rownames(dataRef))
potHetName<- as.numeric(rownames(dataRefDirty))
falsePosIndex <- which(! potHetName %in% trueHetName)
falsePosPosition<- as.integer(potHetName[falsePosIndex])



# cnvTF<- vector(length=length(falsePosPosition))
# for(i in 1:length(falsePosPosition) ){
#     index<- falsePosPosition[i]
#     if( any(which(index >= uniqueCnvRegion[,1] & index <= uniqueCnvRegion[,2]  ))  ){
#         cnvTF[i] <- TRUE
#     }
# }

cnvTF<- sapply(falsePosPosition, function(x){
        if( any((x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  ))  ){
            return(TRUE)
        }
        return(FALSE)
    })

summary(cnvTF)
cnvCount_T<- sum(cnvTF)
cnvCount_F<- sum(!cnvTF)

cat(cnvCount_F,  cnvCount_T, formatC(cnvCount_T/length(cnvTF)), "\n")

# 
trueHetName<- as.numeric(trueHetName)
# # allTF<- vector(length=length(trueHetName))
# # for(i in 1:length(trueHetName) ){
# #     index<- trueHetName[i]
# #     if( any(which(index >= uniqueCnvRegion[,1] & index <= uniqueCnvRegion[,2]  ))  ){
# #         allTF[i] <- TRUE
# #     }
# # }
# 
allTF<- sapply(trueHetName, function(x){
        if( any((x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  ))  ){
            return(TRUE)
        }
        return(FALSE)
    })


summary(allTF)
count_T<- sum(allTF)
count_F<- sum(!allTF)

cat(count_F,  count_T, formatC(count_T/length(allTF)), "\n")
 

combinedCnvCount<- matrix(c(table(cnvTF), table(allTF)), nrow=2, byrow=T)
r<- fisher.test(combinedCnvCount)



cat(paste0(fullTitle, " & FH"), " & ",
    paste(cnvCount_F,  cnvCount_T, formatC(cnvCount_T/length(cnvTF)), sep=" & "),
    " & " ,formatC(r$p.value), " \\\\ " , file=cnvResultFileName, append=TRUE, fill=TRUE)


cat(paste0("" , " & TH"), " & ",
    paste(count_F,  count_T, formatC(count_T/length(allTF)), sep=" & "),
    " & \\\\ \\hline" , file=cnvResultFileName, append=TRUE, fill=TRUE)
# unique( cbind( rep(1:2, 60), rep(1:5, 24)) )


}
cat(sufix, file=cnvResultFileName, fill=T, append=T)



################################################################################
################################################################################
################################################################################

dl<- vector(lengtth(true

dl<- sapply(trueHetName, function(x){
        which(x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  )
#         if( any((x >= uniqueCnvRegion[,1] & x <= uniqueCnvRegion[,2]  ))  ){
#             return(TRUE)
#         }
#         return(FALSE)
    })

dl2<- sapply(dl, function(x){
    dd[x]
})


dl_type<- sapply(dl, function(x){
    return(list(cnvChromosome[x,6]))
})
a=sapply(dl_type,summary)
apply(a,1,sum)


dlx<- sapply(dl2, function(x){
    sum(x>10000)/length(x)
})





##
##
site<- list()
for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)
chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)
setwd(fullPath)

hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE)


trueHetName<- as.numeric(rownames(dataRef))
potHetName<- as.numeric(rownames(dataRefDirty))
falsePosIndex <- which(! potHetName %in% trueHetName)



siteData<- as.data.frame(matrix(nrow=NROW(dataRefDirty), ncol=2))
colnames(siteData)<- c("site", "True/Potential_Het")
siteData[,1]<- potHetName
siteData[,2]<- "T"
siteData[falsePosIndex,2]<- "P"
site[[chromosomeIndex]]<- rbind(site[[chromosomeIndex]], siteData)
# write.table(siteData, file=paste0(pwd,"siteData_",subName), row.names=F)

}


dim(unique(site[["21"]]))
dim((site[["21"]]))

dim(unique(site[["10"]]))
dim((site[["10"]]))


a<-order(unique(site[["10"]])[,1])
site[["10Sort"]] <- unique(site[["10"]])[a,]
write.table(site[["10Sort"]], file=paste0(pwd,"CEU_SiteData_Chromosome10"), row.names=F)

a<-order(unique(site[["21"]])[,1])
site[["21Sort"]] <- unique(site[["21"]])[a,]
write.table(site[["21Sort"]], file=paste0(pwd,"CEU_SiteData_Chromosome21"), row.names=F)




