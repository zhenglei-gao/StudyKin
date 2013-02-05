# Trial           : guitest
# File name       : guitest IRLS SFO  parent.r
# Target path     : D:\Projects\Kinetics\KINGUI\demo\guitest
# Created         : on 12 Aug 2011
#                   at 03:22
#                   by gbbfx on ADEMONC5902(2CPUs)
# KinGUII version : 2.2011.811.11535
# Comments        :

#
# Settings
#

options(warn=-1)
rm(list=ls())
getwd()

sourceDir <- function(path, trace = TRUE, loadSourcefun=TRUE,loadRdata=FALSE,...) {
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
# 
sourceDir('C:/Users/gbbfx/Work/Projects/Kinetics/KINGUI/package/workingcopy/Sourcefunctions/')
#
guitest <- mkinmod.full(
  parent = list(
       time = c(     0,      3,      7,     14,     30,     62,     90,    118),
    residue = c(101.24,  99.27,  90.11,  72.19,  29.71,   5.98,   1.54,  NA),
     weight = c(     1,      1,      1,      1,      1,      1,      1,      1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.040,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 100.15,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf)),
          inpartri='default',outpartri='default' )

#
# Fit and optimizer
#

  Fit    <- IRLSkinfit.full(
            guitest,
               plot      = TRUE,
               quiet     = TRUE,
               ctr       = kingui.control(
                               method = 'solnp',
                            submethod = 'Port',
                              maxIter = 100,
                            tolerance = 1E-06,
                            odesolver = 'lsoda'),
            irls.control = list(
                              maxIter = 10,
                            tolerance = 0.001))

#
# Output
#

list2ascii(
            summary(Fit, cov = TRUE,version=version),
            'Focus_2006_A.kgo')

kingraph(Fit,
            'Focus_2006_A.kgg')

kinplot(Fit)
#
# End of R-script
#

