##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param lib
##' @param pkg
##' @return
##' @author
##'
##version <- read.table('versioninfo.txt',header=TRUE,sep='\t')[1,1]
.First.lib <- function(lib, pkg){

 if (!require(mkin, quietly = TRUE))
   warning('mkin could not be loaded')
 if (!require(minqa, quietly = TRUE))
     warning('minqa should be loaded')
 if (!require(Rsolnp, quietly = TRUE))
   warning('Rsolnp be loaded')
}
.First <- function(){

 if (!require(mkin, quietly = TRUE))
   warning('mkin could not be loaded')
 if (!require(minqa, quietly = TRUE))
     warning('mminqa should be loaded')
 if (!require(Rsolnp, quietly = TRUE))
   warning('Rsolnp be loaded')
 if (!require(optimx, quietly = TRUE))
   warning('optimx be loaded')
 options(warn=-1)
}
.First()
