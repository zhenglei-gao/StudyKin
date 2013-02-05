readversion <- function(path,default="1.2011.701.11155", ...) {

    version <- try(as.character(read.table(paste(path,'versioninfo.txt',sep=''),header=TRUE,sep='\t')[1,1]),silent=TRUE)
    if(class(version)=='try-error') version <- default
    return(version)
 }
