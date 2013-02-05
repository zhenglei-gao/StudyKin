# Trial           : complex
# File name       : complex IRLS SFO  parent.r
# Target path     : D:\Projects\Kinetics\KINGUI\demo\guitest\complex
# Created         : on 12 Aug 2011
#                   at 03:35
#                   by gbbfx on ADEMONC5902(2CPUs)
# KinGUII version : 2.2011.811.11535
# Comments        :

#
# Settings
#

options(warn=-1)
rm(list=ls())
library(mkin)
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


##sourceDir('D:/Projects/Kinetics/KINGUI/GUI/KinGUII/R-2.12.1/Sourcefunctions/')
sourceDir('D:/Projects/Kinetics/KINGUI/package/Sourcefunctions/')
version <- readversion('D:/Projects/Kinetics/KINGUI/package/Sourcefunctions/')
complex0 <- mkinmod.gui(
     parent = list(
       time = c(   0,    1,    3,    7,   14,   30,   62,  100),
    residue = c(93.2, 89.4, 79.7, 61.1, 48.2, 15.9,  6.5,    6),
     weight = c(   1,    1,    1,    1,    1,    1,    1,    1),
                      to = c("A1","B1","C1"),
         FF = list(ini   = c(0.1,0.1,0.1),
                   fixed = c(0,0,0),
                   lower = c(0.0,0.0,0.0),
                   upper = c(1.0,1.0,1.0)),
                   sink  = FALSE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf)),
         A1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,    NA,  0.55,  6.87, 17.08, 21.68, 15.77, 13.63),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                      to = c("A2"),
         FF = list(ini   = c(0.1),
                   fixed = c(0),
                   lower = c(0.0),
                   upper = c(1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         A2 = list(
       time = c(   0,    1,    3,    7,   14,   30,   62,  100),
    residue = c(  NA, 0.55, 1.41, 0.55, 1.29, 1.95, 3.54, 3.86),
     weight = c(   1,    1,    1,    1,    1,    1,    1,    1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         B1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,    NA,    NA,  0.55,  2.31, 15.76,  6.36,  3.74),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         C1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,  0.55,   3.2,  5.46, 12.55, 10.45,  4.74,  4.33),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)))
Fit0    <- mkinfit.gui(
            complex0,
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

list2ascii(
            summary(Fit0, cov = TRUE,version=version),
            'schaefer07_complex_case_default.kgo')
#
#
# Residue data, pathway and kinetics
#

 complex <- mkinmod.full(
     parent = list(
       time = c(   0,    1,    3,    7,   14,   30,   62,  100),
    residue = c(93.2, 89.4, 79.7, 61.1, 48.2, 15.9,  6.5,    6),
     weight = c(   1,    1,    1,    1,    1,    1,    1,    1),
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
         A1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,    NA,  0.55,  6.87, 17.08, 21.68, 15.77, 13.63),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                      to = c("A2"),
         FF = list(ini   = c(0.1),
                   fixed = c(0),
                   lower = c(0.0),
                   upper = c(1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         A2 = list(
       time = c(   0,    1,    3,    7,   14,   30,   62,  100),
    residue = c(  NA, 0.55, 1.41, 0.55, 1.29, 1.95, 3.54, 3.86),
     weight = c(   1,    1,    1,    1,    1,    1,    1,    1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         B1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,    NA,    NA,  0.55,  2.31, 15.76,  6.36,  3.74),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
         C1 = list(
       time = c(    0,     1,     3,     7,    14,    30,    62,   100),
    residue = c(   NA,  0.55,   3.2,  5.46, 12.55, 10.45,  4.74,  4.33),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.1,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)),
                    ## inpartri='water-sediment',outpartri='water-sediment' )
inpartri='default',outpartri='water-sediment' )
                   ## inpartri='default',outpartri='default'     )


# Fit and optimizer
#

  Fit    <- mkinfit.full(
            complex,
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
            'schaefer07_complex_case_water_sediment.kgo')

## list2ascii(
##             summary(o1, cov = TRUE,version=version),
##             'schaefer07_complex_case_default.kgo')
## list2ascii(
##             summary(o2, cov = TRUE,version=version),
##             'schaefer07_complex_case_water_sediment.kgo')

## list2ascii(
##             summary(o4, cov = TRUE,version=version),
##             'schaefer07_complex_case_water_sediment_mcmc.kgo')
   kingraph(Fit,
            'schaefer07_complex_case_default.kgg')


#
# End of R-script
#

