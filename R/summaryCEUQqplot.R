#suppressPackageStartupMessages(
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")

dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10

latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"
pwd <- "/home/steven/Postdoc2/Project_MDM/CEU/"
setwd(pwd)

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


plotqq<- function(z, ff, outerText){
    mains = c("Reference Allele", "Alternate Allele", "Error")

    for(i in 1:3) {
        qqplot(z[,i],ff[,i],xlim=c(0,1),ylim=c(0,1),xlab="Estimated Frequency",ylab="Observed Frequency",main=mains[i])
        abline(0,1)
    #   print(ad.stat.k(ff[,i],z[,i]))
    }
    mtext(outerText, side=3, outer=T, line=-2, cex=2)

}


p<- 8
# for (p in 2:4){
# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# setwd(fullPath)
hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE)

rowSumDataRef<- rowSums(dataRef)
colSumDataRef<- colSums(dataRef)
rowSumDataRefDirty<- rowSums(dataRefDirty)

freqDataRef<- dataRef/rowSumDataRef
freqDataRefDirty<- dataRefDirty/rowSumDataRefDirty


maxModel<- extractMaxModel(fullPath)
whichIsDirty <- grepl("_[0-9]D",names(maxModel))
header<- gsub("hets_" , "", names(maxModel) )

qqplotFile<- file.path(latexDir, paste0("qqPlots_", subName, ".pdf") )


pdf(file=qqplotFile, width=12, height=6)
par(mai=c(0.6,0.7,0.2,0.1), mfrow=c(1,3), 
    cex.main=1.2^4,cex.lab=1.2^2, 
    omi=c(0,0,0.5,0) )

# mains = c("Reference Allele", "Alternate Allele", "Error")

#simple multinomial (no bias)
p1<- (colSumDataRef[1]+colSumDataRef[2])/2
prob <- c(p1, p1, colSumDataRef[3]) /sum(colSumDataRef)
ll <- sum(log(prob)*colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, "\n**** Multinomial ****\n")


#multinomial (with ref bias)
prob <- colSumDataRef / sum(colSumDataRef)
ll <- sum(log(prob)* colSumDataRef)
# cat(sprintf("  ll = %0.16g\n", ll));#print(prob)

b <- rmultinomial(length(rowSumDataRef), rowSumDataRef, prob)
z <- b/rowSums(b)
plotqq(z, freqDataRef, "\n**** Biased Multinomial ****\n")

for( m in 1:length(maxModel)) {

    #dirichlet-multinomial mixture
    mModel<- maxModel[[m]]
    if(! whichIsDirty[m]){
        b <- rmdm(length(rowSumDataRef), rowSumDataRef, mModel$f,
                    phi=mModel$params[,1], p=mModel$params[,2:4])
        z <- b/rowSums(b)
        ff<- freqDataRef    
    } else{
#         b <- rmdm(length(rowSumDataRefDirty), rowSumDataRefDirty, mModel$f,
#                     phi=mModel$params[,1], p=mModel$params[,2:4])
#         ff<- freqDataRefDirty
        next
    }
    
    plotqq(z, ff, paste0("\n**** Dirichlet-Multinomial ", header[m] ," ****\n") )

}


dev.off()

# } # match (p in 1:length(subNameList) ){

########################################
# originally run it with dm.R
# b <- rmdm(length(n),n,m3$param.p,m3$param.a)
# #NOTE: difference between $param $parmas
# # mdmParams should be $params but $params don't have class
# class(m3$params)<- "mdmParams" ## Wont work if class not forwords
# b <- rmdm(length(n),n,m3$f,m3$params)
