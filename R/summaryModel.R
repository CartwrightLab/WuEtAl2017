#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- TRUE
suppressPackageStartupMessages( source("./summaryFunctions.R") )
source("./summarySetup.R")

msTable<- vector(mode="list", length=length(subNameList)*2)
AicThIndex <- c(4,3,3,4,4,4,3,4)
BicThIndex <- BICIndexList/2

AicPhIndex <- c(5,6,6,6,6,6,3,6)
BicPhIndex <- c(3,4,4,5,4,5,3,6)


## latex parameter table    
paramSummaryToTex<- function(x, selectCol){
    newOrder<-rev(order(x[,8]))
    x<- x[newOrder,selectCol]
    if(!is.matrix(x)){
        x<- matrix(x,ncol=length(x))
    }
    index<- c(ncol(x)-1,ncol(x))
    x2 = matrix(sprintf("%#0.4g", x),nrow=nrow(x))
    # x2<- cbind( formatC(x[,-index, drop=F], digits=4, width=8, format="f"),
    #             formatC(x[,index, drop=F], digits=1, width=8, format="f") )
    s<- apply(x2,1,function(y){
        paste(y, collapse=" & ")
    })
    return(s)
}

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

    ##Summary
    modelParametersSummary<- sapply(maxModel, function(x){
        x<- cbind(x$params, mdmAlphas(x$params), prop=x[["f"]] )
        x<- x[,c(2:7,1,8)]
        x<- matrix(x, ncol=8)
        colnames(x)<- c("pi_ref", "pi_alt", "pi_err", "alpha_ref", "alpha_alt", "alpha_err", "phi", "prop")
        
        return(x)
    })


    propEM<- sapply(maxModel, function(x){x$f})
    propML<- getMaxLikelihoodProp(maxLikelihoodTable)

    modelParametersSummary2 <- modelParametersSummary
    propResultTable<- vector("list", length(propEM))
    names(propResultTable)<- names(propEM)
    for(i in 1:length(propML)){
        diff <- propEM[[i]] - propML[[i]]
        propResultTable[[i]]<- cbind(propEM[[i]], propML[[i]], diff)
        colnames(propResultTable[[i]])<- c("EM", "ML", "Diff")
        modelParametersSummary2[[i]] <- cbind(modelParametersSummary[[i]], "ML"=propML[[i]])
        
    }


    modelParametersSummary3<- vector(length=length(modelParametersSummary2)+2, mode="list")
    modelParametersSummary3[3:14] <- modelParametersSummary2
    names(modelParametersSummary3)<- c("Multinomial", "Mulitnomial w/bias", names(modelParametersSummary2))

    colSumDataRef<- colSums(dataRef)

    ## simple multinomial (no bias)
    p1<- (colSumDataRef[1]+colSumDataRef[2])/2
    prob_MN <- c(p1, p1, colSumDataRef[3]) / sum(colSumDataRef)
    modelParametersSummary3[[1]]<- matrix(c(prob_MN,1/0,1/0,1/0,0,1,1), nrow=1, byrow=T)

    ## multinomial (with ref bias)
    prob_MNB <- colSumDataRef / sum(colSumDataRef)
    modelParametersSummary3[[2]]<- matrix(c(prob_MNB,1/0,1/0,1/0,0,1,1), nrow=1, byrow=T)
      
    # formatC(x, digits=3, width=8, format="g")
    # sprintf("%8.3g" , x)
    mpsName<- names(modelParametersSummary3)
    latexTable<- sapply(modelParametersSummary3, function(x){
        s<- paramSummaryToTex(x, 1:9)
        latex<- paste(s, collapse=" \\\\ \n & & ")
    })
    
    k<- modelParametersSummary3[[BICIndexList[p]+1]]
    msTable[[p]]<- paramSummaryToTex(k, c(1:3,7:9))
    
    k<- modelParametersSummary3[[BicPhIndex[p]*2+2]]
    msTable[[p+length(subNameList)]]<- paramSummaryToTex(k, c(1:3,7:9))
    
    
    header<- names(latexTable)
    header<- gsub("hets_" , "", header)
    header<- gsub("_" , " M", header)
    header[-(1:2)]<- gsub("$", " & TH", header[-(1:2)])
    header<- gsub("D & TH", " & PH", header)
    header = gsub("CEU\\d\\d M", "", header)
    header[1] = "M & TH"
    header[2] = "BM & TH"

    # prefix<- paste0("\\begin{tabular}{|c|ccc|ccc|c|c|c|} \n",
    # "\\hline \\multicolumn{10}{|c|}{Parameter estimates ", fullTitle, "}\\\\ \\hline \n",
    # "Model & $\\pi_{ref}$ & $\\pi_{alt}$ & $\\pi_{err}$ & $\\alpha_{ref}$ & $\\alpha_{alt}$ & $\\alpha_{err}$ & $\\varphi$ &  $\\%$ & ML-\\% \\\\ \\hline")

    # suffix<- "\\end{tabular}"
    
    prefix = ""
    suffix = ""

    latexTableFull <- latexTable
    fileLatexTable <- file.path(latexTableDir, paste0(subName,"_latexTable.tex") )
    
    cat(prefix, file=fileLatexTable, fill=T)
    for(i in 1:length(latexTable) ){
        latexTableFull[i]<- paste0(header[i], " & ", latexTable[[i]], " \\\\ \\hline")
        cat(latexTableFull[i], file=fileLatexTable, fill=T, append=T)
    }
    cat(suffix, file=fileLatexTable, fill=T, append=T)
    # write.table(latexTableFull, file=file.path(fullPath, "latxTable"), quote=F, row.names=F)
    

    ## maxlikelihood table
    modelLikelihood<- sapply(maxModel, function(x){x$ll})

    fileMaxLikelihoodLatexTable <- file.path(latexTableDir, paste0(subName, "_maxLikelihoodLatexTable.tex") )

    # prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|c|c|} \n",
    #     "\\hline \\multicolumn{7}{|c|}{", fullTitle, " } \\\\ \\hline \n",
    #     "Model & TH lnL & PH lnL & TH AIC & PH AIC & TH BIC & PH BIC \\\\ \\hline")
    # suffix<- "\\hline\n\\end{tabular}"

    prefix = paste0("\\multirow{6}{*}{\\rotatebox[origin=c]{90}{", fullTitle, "}}")
    suffix = ""

    maxLiTable<- matrix(ncol=6, nrow=length(modelLikelihood)/2 )
    numFreeP<- seq(3,by=4,length=6)
    coefBIC<- c(log(NROW(dataRef)), log(NROW(dataRefDirty)) )
    for(i in 1:NROW(maxLiTable) ){
        t_li<- c(modelLikelihood[i*2-1], modelLikelihood[i*2])
        t_AIC<- -2*t_li + 2*numFreeP[i]
        t_BIC<- -2*t_li + coefBIC*numFreeP[i]
        maxLiTable[i,]<- c(t_li, t_AIC, t_BIC)
    }
    minIndex<- apply(maxLiTable,2,which.min)
    maxLiTable<- formatC(maxLiTable,  digits=2, width=12, format="f")
    for(i in 3:NCOL(maxLiTable)){
        maxLiTable[minIndex[i], i]<-paste0("\\textbf{",maxLiTable[minIndex[i], i], "}")
    }

    maxLikelihoodLatex<- apply(maxLiTable,1,function(y){
        paste(y, collapse=" & ")
    })


    cat(prefix, file=fileMaxLikelihoodLatexTable, fill=T)
    for(i in 1:length(maxLikelihoodLatex) ){
        latex<- paste0(" & ", i , " & ", maxLikelihoodLatex[i], " \\\\ ")
        cat(latex, file=fileMaxLikelihoodLatexTable, fill=T, append=T)
    }
    cat(suffix, file=fileMaxLikelihoodLatexTable, fill=T, append=T)


}  ## end for loop for(p in length(subNameList):1 )



################################################################################
##### MS Tables
################################################################################

prefix<- paste0("\\begin{tabular}{|c|c|ccc|c|cc|} \n",
     "\\hline \n",
"Dataset &  BIC & $\\pi_{ref}$ & $\\pi_{alt}$ & $\\pi_{err}$ & $\\varphi$ &  $\\rho$ & ML-$\\rho$ \\\\ \\hline")

suffix<- "\\end{tabular}"


fileLatexMSTable <- file.path(latexTableDir, paste0(projectName,"_maxLikelihoodLatexTable.tex") )

cat(prefix, file=fileLatexMSTable, fill=T)
# for(p in length(subNameList):1 ){
for(p in c(rev(grep("C21", subNameList)), rev(grep("C10", subNameList))) ){
    fullTitle<- fullTitleList[p]
   
    header<- sub(" ", " TH ", fullTitle)
#     formattedText<- paste0( paste(header, AicThIndex[p], BicThIndex[p], sep=" & "),
    formattedText<- paste0( paste(header, BicThIndex[p], sep=" & "),
            " & ", msTable[[p]][1] , " \\\\ \n",
            paste(" &   & ",   msTable[[p]][-1], " \\\\", collapse=" \n"),
            " \\hline")
    cat(formattedText, file=fileLatexMSTable, fill=T, append=T)
    
}
cat("\n%%PH dataset\n", file=fileLatexMSTable, fill=T, append=T)

offset <- length(subNameList)
for(p in c(rev(grep("C21", subNameList)), rev(grep("C10", subNameList))) ){
    fullTitle<- fullTitleList[p]
   
    header<- sub(" ", " PH ", fullTitle)
    formattedText<- paste0( paste(header, BicPhIndex[p], sep=" & "),
            " & ", msTable[[p+offset]][1] , " \\\\ \n",
            paste(" &   & ",   msTable[[p+offset]][-1], " \\\\", collapse=" \n"),
            " \\hline")
    cat(formattedText, file=fileLatexMSTable, fill=T, append=T)
    
}
cat(suffix,  file=fileLatexMSTable, fill=T, append=T)

##############################################################################
##### SNP Count Table
##############################################################################

## SNP count
if(isCEU){
    fileSnpCountLatexTable <- file.path(latexTableDir, paste0("snpCountLatexTable.tex") )
} else {
    fileSnpCountLatexTable <- file.path(latexTableDir, paste0("snpCountLatexTable_CHM1.tex") )
}
# "\tDataset & Individual caller only & Trio caller only  & Both callers (PH) &",
#     "True heterozygotes (TH)  \\\\ \n",
prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|}\n",
    "\t\\hline \n",
    "\tDataset & Individual & Trio & Both &",
    "True   \\\\ \n",
    "\tDataset & caller only & caller only  & callers (PH) &",
    "heterozygotes (TH) \\\\ \n",
    "\t\\hline \n")
suffix<- "\\hline\n\\end{tabular}"


cat(prefix, file=fileSnpCountLatexTable, fill=T)
for(p in length(subNameList):1 ){

    subName<- subNameList[p]
    fullTitle<- fullTitleList[p]
    fullPath <- file.path(dataDir, subName)
   
    temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
    dataFull <- temp$dataFull

    name2<- gsub("_" , " ", subName)
    name2<- gsub("_C" , " Chr", subName)

    countThree<- table(dataFull$callby)[c(1,3,2)]
    countTH<- sum( (dataFull$callby == 2 & (dataFull$snp == 1 & dataFull$snpdif == 0)) ) 
#     percentTH<- formatC(countTH/countThree[3], digits=3, width=8, format="g")
    snpCount <- c(name2, countThree, countTH)#, percentTH)  
 
    snpLatex <- paste(snpCount, collapse=" & ")
    snpLatex <- paste(snpLatex ,"\\\\")
    cat(snpLatex, file=fileSnpCountLatexTable, fill=T, append=T)
}
cat(suffix, file=fileSnpCountLatexTable, fill=T, append=T)
