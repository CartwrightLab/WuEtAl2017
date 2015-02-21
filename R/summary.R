
pid <- Sys.getpid()
args <- commandArgs(trailingOnly=TRUE)
dirty_data <- FALSE

if(length(args) < 3) {
	message("  Usage: search-gt-mdm (num of trials) (num of compoents[d]) (data file).")
	message("  To do a search using dirty data, append 'd' to the number of components.")
	quit(status=1)
}

nn <- as.numeric(args[1])
k <- args[2]
dirty_data <- grepl("^\\d+[Dd]$",k)
if(dirty_data) {
	k <- as.numeric(sub("[Dd]$","",k))
} else {

pwd <- "/home/steven/Postdoc2/Project_MDM/run/"
setwd(pwd)

sub_folder <- c("2010/multiallelic/base_count/",
"2010/original/base_count/",
"2011/multiallelic/base_count/",
"2012/multiallelic/base_count/",
"2012/original/base_count/",
"2013/multiallelic/base_count/",
"2013/original/base_count/")

full_path <- paste(pwd, sub_folder, sep="")

for(p in 1:length(full_path){

	setwd(full_path[p])
	list_rdata<- list.files(path=full_path[p], pattern="RData") 
	max_data<- list()


	for(i in 1:length(list_rdata)){

		load(list_rdata[i])
		data<- get(resName)
		if(length(data) != 0){
			index <- which.max(sapply(data, function(x) {x$ll} ))
			max_data[[i]]<- data[[index]]
		}
		else{
			max_data[[i]]<- NULL
		}
			

		sapply(data, function(x) {x$params} )

	}
	sapply(max_data, function(x){x$params})
	
	
}















