
sourceDir <- function(path, trace = TRUE, loadSourcefun=TRUE,loadRdata=FALSE,...) {
    ## modified from the sourceDir function in the example of source help files.
    if(loadSourcefun==TRUE){
        for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
            if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
        }
    }
   if(loadRdata==TRUE){
       ##for (nm in list.files(paste(path,'data',sep=''), pattern = "*.Rdata")) {
       ##browser()
       for (nm in list.files(path, pattern = "*.Rdata")) {
           if(trace) cat(nm,":")
           ##load(file.path(path, nm), env=parent.frame(),...) ## also works
           load(file.path(path, nm), env=.GlobalEnv,...)
           if(trace) cat("\n")
       }

   }
 }
