## R --vallina filename
scriptname <- sub("^--file=",'',grep("^--file=",commandArgs(),value=TRUE)[1])
if(is.na(scriptname)){
    isRscriptMode<- FALSE
    cat("Local mode", fill=T)
} else {
    isRscriptMode<- TRUE
    cat("Running script ", scriptname, fill=T)
}

.Last <- function(){ 
    if(isRscriptMode){
        cat("Done script ", scriptname, fill=T)
    }
}

# suppressPackageStartupMessages( source("./summaryFunctions.R") )


dirtyData <- FALSE
upperLimit <- 150
lowerLimit <- 10
loadData<- TRUE


latexTableDir<- "../tables/"
figureDir<- "../supfigures"
cnvDir <- "../data/CNV/"
lcDir <- "../data/LowComp/"
gatkDir<- "../data/GATK/"

BICIndexList<- c(2,2,2,2,2,2,2,3)*2


if(isCEU){
    dataDir<- "../data/CEU/"
    subNameList<- c(
        "CEU10_C10", "CEU10_C21",
        "CEU11_C10", "CEU11_C21",
        "CEU12_C10", "CEU12_C21",
        "CEU13_C10", "CEU13_C21"
    )

#     fullTitleList<- c(
#         "CEU 2010 Chromosome 10", "CEU 2010 Chromosome 21",
#         "CEU 2011 Chromosome 10", "CEU 2011 Chromosome 21",
#         "CEU 2012 Chromosome 10", "CEU 2012 Chromosome 21",
#         "CEU 2013 Chromosome 10", "CEU 2013 Chromosome 21"
#     )
    fullTitleList<- c(
        "CEU10 Chr10", "CEU10 Chr21",
        "CEU11 Chr10", "CEU11 Chr21",
        "CEU12 Chr10", "CEU12 Chr21",
        "CEU13 Chr10", "CEU13 Chr21"
    )
    projectName<- "CEU"
    
} else {
    dataDir<- "../data/CHM1/"
    subNameList<- c(
        "CHM1_C10", "CHM1_C21"
    )

#     fullTitleList<- c(
#         "CHM1 Chromosome 10", "CHM1 Chromosome 21"
#     )
    fullTitleList<- c(
        "CHM1 Chr10", "CHM1 Chr21"
    )
    projectName<- "CHM1"
}

 