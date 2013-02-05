sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
 }
sourceDir('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')
version <- readversion('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')
################################
complex <- mkinmod.full(
  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
  A1 = list(type = "SFO", to = "A2"),
  B1 = list(type = "SFO"),
  C1 = list(type = "SFO"),
  A2 = list(type = "SFO"),
  inpartri='default',
  outpartri='default',
  data=schaefer07_complex_case,
  weight=NULL)

####################################
 complex <- mkinmod.full(
parent = list(time = c(0, 1, 3, 7,14,30,62,100),
 residue = c(93.2,89.4,79.7,61.1,48.2,15.9,6.5,6),
 weight = c(   1, 1, 1, 1, 1, 1, 1, 1),
 to = c("A1","B1","C1"),
 FF = list(ini   = c(0.1,0.1,0.1),
           fixed = c(0,0,0),
           lower = c(0.0,0.0,0.0),
           upper = c(1.0,1.0,1.0)),
  sink  = TRUE,
  type = "SFO",
  k = list(ini   = 0.1,
           fixed = 0,
           lower = 0.0,
           upper = Inf),
  M0 = list(ini   = 0.1,
            fixed = 0,
            lower = 0.0,
            upper = Inf)),
 A1 = list( time = c(0,1,3,7, 14, 30, 62,100),
 residue = c(NA, NA,0.55,6.87,17.08,21.68,15.77,13.63),
  weight = c(1,1,1,1,1,1,1,1),
      to = c("A2"),
      FF = list(ini   = c(0.1), fixed = c(0),lower = c(0.0), upper = c(1.0)),
      sink  = TRUE,
      type = "SFO",
      k = list(ini=0.1,fixed=0,lower=0.0,upper=Inf),
     M0 = list(ini= 0,fixed=1,lower=0.0,upper=Inf)),
 A2 = list(time = c(   0, 1, 3, 7,14,30,62,100),
     residue = c(NA,0.55,1.41,0.55,1.29,1.95,3.54,3.86),
     weight = c(   1, 1, 1, 1, 1, 1, 1, 1),
     sink  = TRUE,
     type = "SFO",
     k = list(ini=0.1,fixed=0,lower=0.0,upper=Inf),
     M0 = list(ini=0, fixed=1,lower=0.0,upper=Inf)),
 B1 = list(time = c(    0,1,3,7, 14, 30, 62,100),
    residue = c(NA,NA,NA,0.55,2.31,15.76,6.36,3.74),
    weight = c(1,1,1,1,1,1,1,1),
    sink  = TRUE,
    type = "SFO",
    k = list(ini= 0.1,fixed=0,lower=0.0,upper=Inf),
    M0 = list(ini=0,fixed=1, lower=0.0, upper=Inf)),
 C1 = list(time = c(0,1,3,7, 14, 30, 62,100),
     residue = c(NA,0.55,3.2,5.46,12.55,10.45,4.74,4.33),
     weight = c(1,1,1,1,1,1,1,1),
     sink  = TRUE,
     type = "SFO",
     k = list(ini=0.1,fixed=0, lower=0.0, upper=Inf),
     M0 = list(ini= 0,fixed=1,lower=0.0, upper=Inf)))
###############################################################

a <- mkinmod.full(
       parent = list(type = "SFO", to = c("A1", "B1", "C1"),
                     sink = FALSE),
       A1 = list(type = "SFO", to = "A2"),
       B1 = list(type = "SFO"),
       C1 = list(type = "SFO"),
       A2 = list(type = "SFO"),
                  inpartri='water-sediment',
                  outpartri='water-sediment',
                  data=schaefer07_complex_case,
                  weight=NULL)
###############################################################
 Fit    <- mkinfit.full(
            complex,
               plot      = TRUE,
               quiet     = TRUE,
               ctr
                         = kingui.control(
                               method = 'solnp',
                            submethod = 'Port',
                              maxIter = 100,
                            tolerance = 1E-06,
                            odesolver = 'lsoda'),
            irls.control = list(
                              maxIter = 10,
                            tolerance = 0.001))
kinplot(Fit)
dev.print(pdf,'NLSfit.pdf')
###############################################################


Fit    <- IRLSkinfit.full(
            complex,
               plot      = TRUE,
               quiet     = TRUE,
               ctr       = kingui.control
                            (method = 'solnp',
                            submethod = 'Port',
                              maxIter = 100,
                            tolerance = 1E-06,
                            odesolver = 'lsoda'),
            irls.control = list(
                              maxIter = 10,
                            tolerance = 0.001))
kinplot(Fit)
dev.print(pdf,'IRLSfit.pdf')
#########################################################
list2ascii(summary(Fit, cov = TRUE,version=version),
           'example.kgo')
kingraph(Fit,'example.kgg')
#########################################################

#########################################################

#########################################################
################ Package Generation #####################
#########################################################

setwd('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/dev/doc')
require( highlight )
driver <- HighlightWeaveLatex(box=TRUE)
Sweave( 'kinGUI2.Rnw', driver = driver )
###################################
library(roxygen2)
## example usage
package.skeleton('helloRoxygen',code_files='C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/DevSF/cran/Sourcefunctions/IRLSkinfit.full.R',force=TRUE)
## from formatR package
tidy.source = function(file = choose.files()) {
   exprs = parse(file)
   for (i in 1:length(exprs)) {
       dep = paste(deparse(exprs[i]), collapse = "\n")
       dep = substring(dep, 12, nchar(dep) - 1)
       cat(dep, "\n")
   }
}
### Example of generating calling graphs #####
library(mvbutils)
inner1 <- function() {
  "This is the inner1 function"
}

inner2 <- function() {
  "This is the inner2 function"
}

outer <- function(x) {
   i1 <- inner1()
   i2 <- inner2()
}
foodweb()
help(foodweb)
foodweb(border = TRUE, expand.xbox = 3,
        boxcolor = "#FC6512", textcolor = "black",
        cex = 1.2, lwd=2)
foodweb(where = "package:survival", prune = "survexp",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext("The survexp function foodweb")
########################################################

sourceDir('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')
foodweb(c('IRLSkinfit.full','mkinfit.full','mcmckinfit.full','modFit1','modMCMC','mkinmod.full','completeCompound', "kingui.control"  , "ForwardCalcFF" , "modCost" ,'summary.kingui','summary.mcmckingui','kingraph','chi2err'),ancestors=FALSE,descendents=FALSE, border =FALSE,expand.xbox = 1.2, boxcolor = "grey",textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for package ", bold(KinGUII),sep=' ')))
dev.print(pdf,'callgraph.KinGUII.pdf')
foodweb( c( find.funs(),find.funs("package:FME"),find.funs('package:mkin')), prune = "mkinmod.full",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mkinmod.full),sep=' ')))
dev.print(pdf,'callgraph.mkinmod.full.pdf')
###########################################################
foodweb(c( find.funs(),find.funs("package:FME"),'solnp'),where=1, prune = "IRLSkinfit.full",border = TRUE,expand.xbox = 2, boxcolor = "#FC6512",textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(IRLSkinfit.full),sep=' ')))
dev.print(pdf,'call.IRLSkinfit.pdf')

foodweb(c( find.funs(),find.funs("package:FME"),'solnp'), prune = "mcmckinfit.full",
        border = TRUE,
        expand.xbox = 1.2, boxcolor = "grey",
        textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mcmckinfit.full),sep=' ')))
dev.print(pdf,'call.mcmckinfit.pdf')

foodweb(c( find.funs(),find.funs("package:FME"),'solnp','optimx','nls.lm','nlm','nlminb','bobyqa'), prune = "modFit1",ancestors=FALSE,
        border = FALSE,
        expand.xbox = 1.5, boxcolor = "grey",
        textcolor = "black", cex = 0.9, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(modFit1),sep=' ')))
dev.print(pdf,'call.modFit1.pdf')

foodweb(c( find.funs(),'solnp'), prune = "mkinfit.full",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(mkinfit.full),sep=' ')))
dev.print(pdf,'call.mkinfit.pdf')



foodweb(c('summary.kingui','print.summary.kingui','myformat'), prune = "summary.kingui",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(summary.kingui),sep=' ')))
dev.print(pdf,'call.summary.kingui.pdf')
##
foodweb(c('summary.mcmckingui','print.summary.mcmckingui','myformat' ,find.funs("package:FME"),'mcmc'), prune = "summary.mcmckingui",
        border = TRUE,
        expand.xbox = 1.5, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(summary.mcmckingui),sep=' ')))
dev.print(pdf,'call.summary.mcmckingui.pdf')

##################################################
foodweb(c( find.funs(),find.funs("package:FME"),'solnp'), prune = "modCost",
        border = TRUE,
        expand.xbox = 2, boxcolor = "#FC6512",
        textcolor = "black", cex = 1.0, lwd=2)
mtext(expression(paste("Calling Graph for function ", bold(modCost),sep=' ')))
#########

