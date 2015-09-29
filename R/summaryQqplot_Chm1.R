#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

isCEU <- FALSE
# isCEU <- TRUE


dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10

latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

# if(isCEU){
#     pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
#     subNameList<- c(
#     "CEU10_C10", "CEU10_C21",
#     "CEU11_C10", "CEU11_C21",
#     "CEU12_C10", "CEU12_C21",
#     "CEU13_C10", "CEU13_C21"
#     )
# 
#     fullTitleList<- c(
#     "CEU 2010 Chromosome 10", "CEU 2010 Chromosome 21",
#     "CEU 2011 Chromosome 10", "CEU 2011 Chromosome 21",
#     "CEU 2012 Chromosome 10", "CEU 2012 Chromosome 21",
#     "CEU 2013 Chromosome 10", "CEU 2013 Chromosome 21"
#     )
# } else {
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
# }
setwd(pwd)


p<- 2
p<- 1

# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# # setwd(fullPath)
# hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
# dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU)
# dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

hets_byref<- file.path("base_count_meta_subsample") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData, isCEU=isCEU)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU=isCEU)
n <- rowSums(dataRefDirty)
propRef <- dataRefDirty[,1]/n
oo <- propRef > 0.8
# oo <- lowerLimit <= n & n <= upperLimit 
dataRefProp <- dataRefDirty[oo,]

plotData<- dataRefDirty[,-2]
plotData<- dataRefProp[,-2]
rowSumDataRef<- rowSums(plotData)
colSumDataRef<- colSums(plotData)
# rowSumDataRefDirty<- rowSums(dataRefDirty)

freqDataRef<- plotData/rowSumDataRef
# freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty


maxModel<- extractMaxModel(fullPath, isCEU=isCEU)
whichIsDirty <- grepl("_[0-9]D",names(maxModel))
whichIsP <- grepl("_[0-9]P",names(maxModel))
header<- gsub("_" , "", names(maxModel) )
header<- gsub("D" , "F", header )
header<- gsub("P" , "R", header )

plotTitle <- paste0("qqPlots_", subName, ".pdf") 
qqplotFile<- file.path(latexDir, plotTitle)


pdf(file=qqplotFile, width=12, height=6, title=plotTitle)
par(mai=c(0.8,0.8,0.2,0.1), mfrow=c(1,2), 
    cex.main=1.2^2,cex.lab=1.2, 
    omi=c(0,0,0.5,0) )

# mains = c("Reference Allele", "Alternate Allele", "Error")

# #simple multinomial (no bias)
prob <- colSumDataRef /sum(colSumDataRef)
ll <- sum(log(prob)*colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, "\n**** Multinomial ****\n") #, xlim=c(0.8,1), ylim=c(0.8,1))


#multinomial (with ref bias)
# prob <- colSumDataRef / sum(colSumDataRef)
# ll <- sum(log(prob)* colSumDataRef)
# # cat(sprintf("  ll = %0.16g\n", ll));#print(prob)
# 
# b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
# z <- b/rowSums(b)
# plotqq(z, freqDataRef, "\n**** Biased Multinomial ****\n")

for( m in 1:length(maxModel)) {

    #dirichlet-multinomial mixture
    mModel<- maxModel[[m]]
    if( whichIsP[m]){
        b <- rmdm(length(rowSumDataRef), rowSumDataRef, mModel$f,
                    phi=mModel$params[,1], p=mModel$params[,c(2,4)])
        z <- b/rowSums(b)
        ff<- freqDataRef    
    } else{
#         b <- rmdm(length(rowSumDataRefDirty), rowSumDataRefDirty, mModel$f,
#                     phi=mModel$params[,1], p=mModel$params[,2:4])
#         ff<- freqDataRefDirty
        next
    }

    plotqq(z, ff, paste0("\n**** Dirichlet-Multinomial ", header[m] ," ****\n") )
#     plotqq(z, ff, paste0("\n**** Dirichlet-Multinomial ", header[m] ," ****\n"), xlim=c(0.8,1), ylim=c(0.8,1) )

}


dev.off()
embedFonts(qqplotFile, options="-DPDFSETTINGS=/prepress")

# } # match (p in 1:length(subNameList) ){

########################################
# originally run it with dm.R
# b <- rmdm(length(n),n,m3$param.p,m3$param.a)
# #NOTE: difference between $param $parmas
# # mdmParams should be $params but $params don't have class
# class(m3$params)<- "mdmParams" ## Wont work if class not forwords
# b <- rmdm(length(n),n,m3$f,m3$params)
