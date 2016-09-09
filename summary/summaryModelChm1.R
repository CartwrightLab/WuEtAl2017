#!/usr/bin/Rscript
## Usage: Rscript --vallina filename.R
isCEU <- FALSE
suppressPackageStartupMessages( source("./summaryFunctions.R") )
source("./summarySetup.R")

dirtyData <- TRUE

msTable<- vector(mode="list", length=length(subNameList)*2)
AicRdIndex <- c(2,2)
AicFdIndex <- c(3,4)
BicIndex <- c(2,2,2,2)



## latex parameter table    
paramSummaryToTex<- function(x, selectCol){
    newOrder<-rev(order(x[,6]))
    x<- x[newOrder,selectCol]
    if(!is.matrix(x)){
        x<- matrix(x,ncol=length(x))
    }
    index<- c(ncol(x)-1,ncol(x))
    x2<- cbind( formatC(x[,-index, drop=F], digits=3, width=5, format="g",flag="#",drop0trailing=FALSE),
                formatC(x[,index, drop=F], digits=3, width=5, format="g",flag="#",drop0trailing=FALSE) )
    s<- apply(x2,1,function(y){
        paste(y, collapse=" & ")
    })
    return(s)
}
# 
# ## latex parameter table    
# paramSummaryToTex<- function(x, selectCol){
#     newOrder<-rev(order(x[,8]))
#     x<- x[newOrder,selectCol]
#     if(!is.matrix(x)){
#         x<- matrix(x,ncol=length(x))
#     }
#     index<- c(ncol(x)-1,ncol(x))
#     x2<- cbind( formatC(x[,-index, drop=F], digits=4, width=8, format="f"),
#                 formatC(x[,index, drop=F]*100, digits=1, width=8, format="f") )
#     s<- apply(x2,1,function(y){
#         paste(y, collapse=" & ")
#     })
#     return(s)
# }


for(p in length(subNameList):1 ){

    subName<- subNameList[p]
    fullTitle<- fullTitleList[p]
    fullPath <- file.path(dataDir, subName)
   
    chromosomeIndex<- gsub(".*_C([0-9]+)", "\\1", subName)

    temp<- loadRawData(fullPath, isCEU, lowerLimit, upperLimit, dirtyData)
    dataFull <- temp$dataFull
#     dataRef <- temp$dataRef
    dataRefDirty <- temp$dataRefDirty

    maxModel<- loadMaxModel(fullPath, subName, loadData, isCEU, isRscriptMode)
    maxLikelihoodTable <- loadMaxLikelihoodTable(fullPath, subName, loadData, maxModel, dataFull, lowerLimit, upperLimit, isCEU, isRscriptMode)
    
    ## filter
    n <- rowSums(dataRefDirty)
    propRef <- dataRefDirty[,1]/n
    oo <- propRef > 0.8
    dataRefProp <- dataRefDirty[oo,]

    ## summary
    modelParametersSummary<- sapply(maxModel, function(x){
        x<- cbind(x$params, mdmAlphas(x$params), prop=x[["f"]] )
        x<- x[,c(2,4,5,7,1,8)]
        x<- matrix(x, ncol=6)
        colnames(x)<- c("pi_ref", "pi_err", "alpha_ref", "alpha_err", "phi", "prop")
        
        return(x)
    })
    modelParametersSummary<- modelParametersSummary[grepl("_[0-9]D|P",names(modelParametersSummary))]


    propEM<- sapply(maxModel, function(x){ x$f })
    propEM<- propEM[grepl("_[0-9]D|P",names(propEM))]
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


    modelParametersSummary3<- vector(length=length(modelParametersSummary2)+1, mode="list")
    modelParametersSummary3[2:13] <- modelParametersSummary2
    names(modelParametersSummary3)<- c("Multinomial", names(modelParametersSummary2))

    

    ## simple multinomial (no bias)
    colSumDataRef<- colSums(dataRefProp)[-2]
    prob_MN <- colSumDataRef /sum(colSumDataRef)
    modelParametersSummary3[[1]] <- matrix(c(prob_MN,1/0,1/0,0,1,1), nrow=1, byrow=T)


    latexTable<- sapply(modelParametersSummary3, function(x){
#         newOrder<-rev(order(x[,6]))
#         x<- x[newOrder,]
#         
#         x<- formatC(x, digits=3, width=8, format="g")
#         x<- matrix(x,ncol=7)
#         s<- apply(x,1,function(y){
#             paste(y, collapse=" & ")
#         })
        s<- paramSummaryToTex(x, 1:7)
        latex<- paste(s, collapse=" \\\\ \n & & ")
    })


    k<- modelParametersSummary3[[BicIndex[p]*2]]
    msTable[[p]]<- paramSummaryToTex(k, c(1:2,5:7))
    names(msTable)[p]<- paste(fullTitle, names(latexTable)[BicIndex[p]*2])
    k<- modelParametersSummary3[[BicIndex[p]*2+1]]
    msTable[[p+length(subNameList)]]<- paramSummaryToTex(k, c(1:2,5:7))
    names(msTable)[p+length(subNameList)]<- paste(fullTitle, names(latexTable)[BicIndex[p]*2+1])

    
    header<- names(latexTable)
    #header<- paste("CHM1", header)
    header<- gsub("_" , "", header)
    header<- gsub("D" , " & FD", header)
    header<- gsub("P" , " & RD", header)
    header[1] = "M & RD"

    # prefix<- paste0("\\begin{tabular}{|c|cc|cc|c|c|c|} \n",
    #     "\\hline \\multicolumn{8}{|c|}{Parameter estimates ", fullTitle, "}\\\\ \\hline \n",
    #     "Model & $\\pi_{ref}$ & $\\pi_{err}$ & $\\alpha_{ref}$ & $\\alpha_{err}$ & $\\varphi$ &  $\\%p$ & ML-\\% \\\\ \\hline")

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
       

    ## maxlikelihood table
    modelLikelihood<- sapply(maxModel, function(x){x$ll})
    modelLikelihood<- modelLikelihood[grepl("_[0-9]D|P",names(modelLikelihood))]


    fileMaxLikelihoodLatexTable <- file.path(latexTableDir, paste0(subName, "_maxLikelihoodLatexTable.tex") )

    # prefix<- paste0("\\begin{tabular}{|c|c|c|c|c|c|c|} \n",
    #     "\\hline \\multicolumn{7}{|c|}{",fullTitle," } \\\\ \\hline \n",
    #     "Model & FD lnL & RD lnL & FD AIC & RD AIC & FD AIC& RD BIC\\\\ \\hline")
    # suffix<- "\\hline\n\\end{tabular}"
   
    prefix = paste0("\\multirow{6}{*}{\\rotatebox[origin=c]{90}{", fullTitle, "}}")
    suffix = ""

    maxLiTable<- matrix(ncol=6, nrow=length(modelLikelihood)/2 )
    numFreeP<- seq(3,by=4,length=6)
    coefBIC<- c(log(NROW(dataRefDirty)), log(NROW(dataRefProp)) )
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

prefix<- paste0("\\begin{tabular}{|c|c|cc|c|cc|} \n",
     "\\hline \n",
"Dataset & BIC & $\\pi_{ref}$  & $\\pi_{err}$ & $\\varphi$ &  $\\rho$ & ML-$\\rho$ \\\\ \\hline")
# "Dataset & AIC & BIC & $\\pi_{ref}$  & $\\pi_{err}$ & $\\varphi$ &  $p$ & ML-$p$ \\\\ \\hline")
suffix<- "\\end{tabular}"


fileLatexMSTable <- file.path(latexTableDir, paste0(projectName,"_maxLikelihoodLatexTable.tex") )

cat(prefix, file=fileLatexMSTable, fill=T)
offset <- length(subNameList)
for(p in length(subNameList):1 ){
    fullTitle<- fullTitleList[p]
   
    header<- sub("CHM1", "RD", fullTitle)
    si<- offset+p
    formattedText<- paste0( paste(header, BicIndex[si], sep=" & "),
            " & ", msTable[[si]][1] , " \\\\ \n",
            paste(" &  & ",   msTable[[si]][-1], " \\\\", collapse=" \n"),
            " \\hline")
    cat(formattedText, file=fileLatexMSTable, fill=T, append=T)
    
}
cat("\n%%FD dataset\n", file=fileLatexMSTable, fill=T, append=T)


for(p in length(subNameList):1 ){

    fullTitle<- fullTitleList[p]
   
    header<- sub("CHM1", "FD", fullTitle)
    si<- p
    formattedText<- paste0( paste(header, BicIndex[si], sep=" & "),
            " & ", msTable[[si]][1] , " \\\\ \n",
            paste(" &   & ",   msTable[[si]][-1], " \\\\", collapse=" \n"),
            " \\hline")
    cat(formattedText, file=fileLatexMSTable, fill=T, append=T)
    
}
cat(suffix,  file=fileLatexMSTable, fill=T, append=T)

