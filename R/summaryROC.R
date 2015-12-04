require(ROCR)
require(OptimalCutpoints)
source("/home/steven/Postdoc2/Project_MDM/MiDiMu/R/summaryFunctions.R")


isCEU <- FALSE
isCEU <- TRUE


dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10

latexDir<- "/home/steven/Postdoc2/Project_MDM/DiriMulti/"

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
} else {
    pwd <- "/home/steven/Postdoc2/Project_MDM/CHM1/"
    subNameList<- c(
    "CHM1_C10", "CHM1_C21"
    )

    fullTitleList<- c(
    "CHM1 Chromosome 10", "CHM1 Chromosome 21"
    )
}
setwd(pwd)


likelihoodToProportion<- function(ml, p){
    eml<- t(t(exp(ml))*p)
    peml<- eml/rowSums(eml)
    return (peml)
}


allAUC<- matrix(nrow=length(subNameList), ncol=6)

p<- 8
# for(p in 1:length(subNameList) ){

subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# setwd(fullPath)
maxModel<- extractMaxModel(fullPath)

hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTableFull.RData")
if ( file.exists(fileMaxLikelihoodTabel) ){
    load(fileMaxLikelihoodTabel)
} else{
    stop(paste0("File does not exist: ", fileMaxLikelihoodTabel))
}

#  refIndex<- (dataFull$pos %in% as.numeric( rownames(dataRef)) )
# snpTfClean <- (dataFull$snp == 1 & dataFull$snpdif == 0)[refIndex]
refDirtyIndex<- (dataFull$pos %in% as.numeric( rownames(dataRefDirty)) )
snpTf <- (dataFull$snp == 1 & dataFull$snpdif == 0)[refDirtyIndex]

whichIsDirty <- grepl("_[0-9]D",names(maxLikelihoodTable))
header<- gsub("hets_CEU13" , "", names(maxLikelihoodTable) )
header<- gsub("D", " components", header)
header<- gsub("_", " ", header)



rocPlotFile<- file.path(latexDir, paste0("rocPlots_", subName, ".pdf") )
pdf(file=rocPlotFile, width=12, height=8, title=rocPlotFile)
par(mar=c(3,3,2,1), mgp=c(1.75, 0.6, 0), #mfrow=c(2,3), 
    cex.main=1.2^3, cex.lab=1.2^2, oma=c(0,0,2.5,0) )

# for( m in 3:length(maxLikelihoodTable)) {
colIndex <- 1
whichDirtyIndex<-which(whichIsDirty)[-1]
legendAUC<- vector()
for(m in whichDirtyIndex){

    maxIndex<- which.max(maxModel[[m]]$f)
    
    ml<- maxLikelihoodTable[[m]]
    classProp<- likelihoodToProportion(ml, maxModel[[m]]$f)
    pred<- prediction(classProp[,maxIndex], snpTf)
    
#     cleanIndex<- refIndex[refDirtyIndex]
#     predClean<- prediction(classProp[cleanIndex,maxIndex], snpTfClean)
#     table(apply(classProp[,],1,which.max))
#     table(apply(classProp[cleanIndex,],1,which.max))

#         ll<-mdm.ll(dataRefDirty, maxModel[[m]]$f, mdmAlphas(maxModel[[m]]$params))
#         pred<- prediction(ll$p.row[,maxIndex], snpTf)

    perfRoc <- performance(pred, "tpr", "fpr")
    perfAuc <- performance(pred,"auc")
    allAUC[p, m/2] <- perfAuc@y.values[[1]]
    auc<- formatC(allAUC[p, m/2])
    legendAUC[m]<- auc
    plot(perfRoc, main="", xlab="1 - specificity", ylab="Sensitivity", col=colIndex, lty=1, lwd=2)
    colIndex <- colIndex+1
    par(new=T)
}

legend(0.7,0.5, legend=paste(header[whichDirtyIndex], legendAUC[whichDirtyIndex]), col=1:5, lwd=3, lty=1, title="Area under ROC curve")
# paste("ROC curve", header[m], "AUC:",auc)

mtext(fullTitle, outer=T, cex=2, line=0)

dev.off()
embedFonts(rocPlotFile, options="-DPDFSETTINGS=/prepress")

# }  ## end for loop 


################################################################################
##### AUC table #####
################################################################################
# header<- names(latexTable)
# header<- gsub("hets_" , "", header)
# header<- gsub("_" , "M", header)
if(NCOL(allAUC) == 6){
    allAUC <- allAUC[,2:6]
}

prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|c|}
    \\hline \\multicolumn{6}{|c|}{Area under ROC curve}\\\\ \\hline
    Dataset & 2 Components & 3 Components & 4 Components & 5 Components & 6 Components  \\\\ \\hline")

sufix<- "\\end{tabular}"

allAucString<- formatC(allAUC, format="f", digits=3)
allAucLatex<- apply(allAucString,1,function(x){paste0(x, collapse=" & ")})

if(isCEU){
    fileAucTable<- file.path(latexDir, paste0("AUC_summary.tex") )
} else{
    fileAucTable<- file.path(latexDir, paste0("AUC_summary_CHM1.tex") )
}
dataName<- gsub("_C", " Chr", subNameList)
cat(prefix, file=fileAucTable, fill=T)
for(i in length(subNameList):1 ){
    
    temp<- paste0(dataName[i], " & ", allAucLatex[i], " \\\\ \\hline")
    cat(temp, file=fileAucTable, fill=T, append=T)
}
cat(sufix, file=fileAucTable, fill=T, append=T)
# write.table(latexTableFull, file=file.path(fullPath, "latxTable"), quote=F, row.names=F)
    

################################################################################
##### Optimal cut for CEU2013 Chromosome 21
################################################################################

p<- 8
subName<- subNameList[p]
fullTitle<- fullTitleList[p]
subFolders <- paste0(subName, "/original/base_count/")
fullPath <- file.path(pwd, subFolders)

# setwd(fullPath)
maxModel<- extractMaxModel(fullPath)

hets_byref<- list.files(path=fullPath, pattern="hets.+byref") 
dataFull <- read.delim(paste(fullPath, hets_byref, sep=""), header=TRUE)
# dataRef<- parseData(dataFull, lowerLimit, upperLimit, dirtyData)
dataRefDirty<- parseData(dataFull, lowerLimit, upperLimit, dirtyData=TRUE, isCEU)

fileMaxLikelihoodTabel <- file.path(fullPath, "maxLikelihoodTableFull.RData")
if ( file.exists(fileMaxLikelihoodTabel) ){
    load(fileMaxLikelihoodTabel)
} else{
    stop(paste0("File does not exist: ", fileMaxLikelihoodTabel))
}

#  refIndex<- (dataFull$pos %in% as.numeric( rownames(dataRef)) )
# snpTfClean <- (dataFull$snp == 1 & dataFull$snpdif == 0)[refIndex]
refDirtyIndex<- (dataFull$pos %in% as.numeric( rownames(dataRefDirty)) )
snpTf <- (dataFull$snp == 1 & dataFull$snpdif == 0)[refDirtyIndex]

whichIsDirty <- grepl("_[0-9]D",names(maxLikelihoodTable))
header<- gsub("hets_CEU13" , "", names(maxLikelihoodTable) )
header<- gsub("D", " components", header)
header<- gsub("_", " ", header)


m<-which(whichIsDirty)[2]
legendAUC<- vector()

maxIndex<- which.max(maxModel[[m]]$f)

ml<- maxLikelihoodTable[[m]]
classProp<- likelihoodToProportion(ml, maxModel[[m]]$f)
pred<- prediction(classProp[,maxIndex], snpTf)

# perfRoc <- performance(pred, "tpr", "fpr")
# perfAuc <- performance(pred,"auc")
# allAUC[p, m/2] <- perfAuc@y.values[[1]]
# auc<- formatC(allAUC[p, m/2])



# # # list all possible methods
# listMethods<- c(
# "MCT", "CB", 
# "MaxSpSe", "MaxProdSpSe", 
# "ROC01", 
# "SpEqualSe", 
# "Youden", "MaxEfficiency", "Minimax", 
# "MaxDOR", "MaxKappa", "PROC01", 
# "NPVEqualPPV", "MaxNPVPPV", 
# "MaxSumNPVPPV", "MaxProdNPVPPV",
# "MinPvalue", "PrevalenceMatching"
# )

        
xx<- as.data.frame(cbind(prop=classProp[,maxIndex], snpTf=snpTf))
listMethods<- c("CB", "ROC01", "Youden")

oo<- optimal.cutpoints("prop", status="snpTf", tag.healthy=0, data=xx, method=listMethods)

# par(mfrow=c(3,3),
# plot(oo, which=1, ylim=c(0,1.1))
# 
# par(mfrow=c(3,3))
# plot(oo, which=1:3)
# ocSummary<- summary(oo)
# ocSummary$p.table
# ocSummary$p.table$Global$Youden[[1]][1:3]


rocPlotFile<- file.path(latexDir, paste0("rocPlotsCutoff_", subName, ".pdf") )
pdf(file=rocPlotFile, width=12, height=8, title=rocPlotFile)


par(mgp=c(1.75, 0.6, 0), #mfrow=c(2,3), 
    cex.main=1.2^3, cex.lab=1.2^2, oma=c(0,0,0,0) )

plot(1 - m[["measures.acc"]][["Sp"]][, 1], m[["measures.acc"]][["Se"]][, 1], 
    xlab = "1-Specificity", ylab = "Sensitivity", 
    main = fullTitle, type = "l", cex.lab = 1.3, cex.axis = 1.3,
    ylim=c(0,1.02) )
# mtext(fullTitle, outer=T, cex=2, line=0)

for(j in 1:length(listMethods)){
m <- oo[[j]][[1]]
cutoff <- m[["optimal.cutoff"]][["cutoff"]][[1]]
sp <- m[["optimal.cutoff"]][["Sp"]][[1]]
x<- 1 - m[["optimal.cutoff"]][["Sp"]][[1]]
se <- m[["optimal.cutoff"]][["Se"]][[1]]
y<- se

points(x, y, pch = 16, cex = 1, col=j)
legend.text <- paste(listMethods[j], " - Cutoff: ", round(cutoff, 3), "\n",
        "Sensitivity: ", round(se, 3), "\n", 
#         paste(rep(" ", nchar(listMethods[j])*2), collapse=""),
        "Specificity: ", round(sp, 3), sep = "")
legend(x, y, col=j, listMethods[j], bty = "n", xjust = 0.05, yjust = 0.5)
legend(0.8, 1-0.2-0.15*j, pch=16, col=j, legend.text, bty = "n")
}

dev.off()
embedFonts(rocPlotFile, options="-DPDFSETTINGS=/prepress")


################################################################################
##### make sure the methods are the same
################################################################################
require("ROCR")

m2<- maxModel[[6]]
x<- dataRefDirty
source("/home/steven/Postdoc2/Project_MDM/RachelCode/dm.R")
ll<-mdm.ll(x,m2$f, mdmAlphas(m2$params))

#
ml<- maxLikelihoodTable[[6]]
pp<- m2$f
#
pr<- ml
p<- pp
prm <- apply(pr,1,max)
prm <- 0
pr <- pr-prm
pr <- t(t(exp(pr))*p)
rs <- rowSums(pr)
rss <- log(rs)+prm
all.equal(sum(rss), ll$ll)
all.equal(rss, ll$ll.row)
all.equal(pr/rs, ll$p.row)

# eml<- t(apply(ml,1,function(x){ exp(x)*p }))
eml<- t(t(exp(ml))*p)
peml<- eml/rowSums(eml)
all.equal(peml, ll$p.row)

classProp<- likelihoodToProportion(ml, p)
all.equal(classProp, ll$p.row)

#     ll<-mdm.ll(dataRefDirty, m2$f, mdmAlphas(m2$params))
#     pred<- prediction(ll$p.row[,3], classicifation)
#     pred<- prediction(peml[,3], snpTf)
        # ll<-mdm.ll(dataRefDirty, m2$f, mdmAlphas(m2$params))
        # pred<- prediction(ll$p.row[,3], classicifation)

        
