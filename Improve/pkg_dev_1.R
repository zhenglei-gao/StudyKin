if ("devtools" %in% loadedNamespaces()) {
  
  stop("You must restart R before installing devtools")
  
}



# -------------------------------------------------------
url <- "https://gist.github.com/raw/4506250/devtools.zip"

temp <- file.path(tempdir(), "devtools.zip")



setInternet2(TRUE)

suppressWarnings(download.file(url, temp, mode = "wb"))

install.packages(temp, repos = NULL)

file.remove(temp)
# -------------------------------------------------------

library(devtools)
has_devel()
#### 
install_github("pryr")
# -------------------------------------------------------
library(roxygen2)
# upper directory of the package
setwd('C:/Users/z.gao/Documents/GitHub/')  
setwd("C:/Projects2013/KinEvalGit/")
setwd("E:/KinEvalGit/")
roxygenize('KineticEval')
check('KineticEval')
document('KineticEval')
build('KineticEval')
install_github("KineticEval","zhenglei-gao")
# --------------------------------------------------------
(.packages())
detach("package:KineticEval")
search()
# --------------------------------------------------------
## Install the package  R CMD INSTALL pkg
## check on cran R CMD check pkg


# -------------------------------------------------------
## knowledges about programming I have learned:
## Example:
## ## In the <environment: R_GlobalEnv>
a <- 1
b <- 2
f1 <- function(){
  c <- 3
  d <- 4
  f2 <- function(P){
    assign("calls", calls+1, inherits=TRUE)
    print(calls)
    return(P+c+d)
  }
  calls <- 0
  v <- vector()
  for(i in 1:10){
    v[i] <- f2(P=0)
    c <- c+1
    d <- d+1
  }
  return(v)
}
f1()
#################
f4 <- function(P){
  assign("calls", calls+1, inherits=TRUE)
  print(calls)
  return(P+c+d)
}
f4_other <- function(P){
  return(P+c+d)
}
f5 <- function(P,liste){
  with(liste,{
    assign("calls", calls+1, inherits=TRUE,envir = sys.frame(-1))
    ##assign("calls", calls+1, inherits=TRUE,envir = sys.frame(-1))
    ##calls <<- calls + 1
    print(calls)
    return(P+c+d)
  }
  )
}

f6 <- function(P,liste){
  attach(liste,pos=2)
  ## assign("calls", calls+1, inherits=FALSE)
    ##assign("calls", calls+1, inherits=TRUE,envir = sys.frame(-1))
  calls <<- calls + 1
  print(paste('HI',i,'HI',calls))
  return(P+c+d)
  
}

##########
f3 <- function(){
  c <- 3
  d <- 4
  calls <- 0
  v <- vector()
  for(i in 1:10){
    ##browser()
    ##v[i] <- f4(P=0) ## or replace here with f5(P=0)
    v[i] <- f5(P=0,liste=as.list(environment()))
    c <- c+1
    d <- d+1
  }
  return(v)
}
######
f3()
##--------------------------------
f7 <- function(P,calls,liste){
  ##calls <<- calls+1
  ##browser()
  assign("calls", calls+1, inherits=TRUE,envir = sys.frame(-1))
  print(calls)
  with(liste,{
    print(paste('with the listed envrionment, calls=',calls))
    return(P+c+d)
  }
  )
}
########
##################
f8 <- function(){
  c <- 3
  d <- 4
  calls <- 0
  v <- vector()
  for(i in 1:10){
    ##browser()
    ##v[i] <- f4(P=0) ## or replace here with f5(P=0)
    v[i] <- f7(P=0,calls,liste=as.list(environment()))
    c <- c+1
    d <- d+1
  }
  f7(P=0,calls,liste=as.list(environment()))
  print(paste('final call number',calls))
  return(v)
}
f8()
calls <- 8
c <- 15
d <- 15
f7(P=0,calls,liste=as.list(environment()))
f7(P=0,calls=1,liste=as.list(environment()))

####################
f2a <- function(P, env = parent.frame()) {
  env$calls <- env$calls + 1
  print(env$calls)
  return(P + env$c + env$d)
}

a <- 1
b <- 2
# same as f1 except f2 removed and call to f2 replaced with call to f2a
f1a <- function(){
  c <- 3
  d <- 4
  calls <- 0
  v <- vector()
  for(i in 1:10){
    v[i] <- f2a(P=0)
    c <- c+1
    d <- d+1
  }
  return(v)
}
f1a()
calls

f2b <- function(P) {
  calls <<- calls + 1
  print(calls)
  a <<- 4
  b <<-9
  return(P + c + d)
}

a <- 1
b <- 2
# same as f1 except f2 removed, call to f2 replaced with call to f2b
#  and line marked ## at the beginning is new
f1b <- function(){
  environment(f2b) <- environment() ##
  c <- 3
  d <- 4
  calls <- 0
  v <- vector()
  for(i in 1:10){
    v[i] <- f2b(P=0)
    c <- c+1
    d <- d+1
  }
  print(paste('now a=', a))
  return(v)
}
f1b()
a
##########


