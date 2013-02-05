# Trial           : Cereals6118Fixed
# File name       : Cereals6118Fixed IRLS SFO  PTZ.r
# Target path     : C:\Emod\KinGUII\WorkingDirectory\Voles\Cereals
# Created         : on 18 Dec 2012
#                   at 08:51
#                   by pfcaa on ADEMONC7403(4CPUs)
# KinGUII version : 2.2012.1030.1351
# Comments        : 

#
# Settings
#

options(warn=-1)
rm(list=ls())
library(mkin)
sourceDir <- function(path, trace = TRUE, ...) {
   for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
 }


sourceDir('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')
version <- readversion('S:/Software & EDV/Models/FOCUS Kinetics/KinGUII/R-2.12.1/SourceFunctions/')


#
# Residue data, pathway and kinetics
#

 Cereals6118Fixed <- mkinmod.full(
        PTZ = list(
       time = c(    0,     0,     1,     1,     3,     3,     7,     7,    16,    16,    21,    21,    28,    28),
    residue = c(1.676, 2.095, 1.362, 2.444, 0.628, 0.506, 0.310, 0.235, 0.075, 0.088, 0.030, 0.033, 0.023, 0.015),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1),
                      to = c("Desthio"),
         FF = list(ini   = c(0.50),
                   fixed = c(0),
                   lower = c(0.0),
                   upper = c(1.0)),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.4000,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 3,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf)),
    Desthio = list(
       time = c(    0,     0,     1,     1,     3,     3,     7,     7,    16,    16,    21,    21,    28,    28),
    residue = c(0.192, 0.219, 0.461, 0.444, 0.375, 0.294, 0.170, 0.175, 0.034, 0.051, 0.020, 0.026, 0.014, 0.011),
     weight = c(    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1),
                   sink  = TRUE,
       type = "SFO",
          k = list(ini   = 0.8000,
                   fixed = 0,
                   lower = 0.0,
                   upper = Inf),
         M0 = list(ini   = 0,
                   fixed = 1,
                   lower = 0.0,
                   upper = Inf)))

#
# Fit and optimizer
#

sourceDir('C:/Projects/KinGui2/workingcopy_2012_June/Sourcefunctions/')

  Fit    <- IRLSkinfit.full(
            Cereals6118Fixed,
               plot      = TRUE,
               quiet     = TRUE,
               ctr       = kingui.control(
                               method = 'Marq',
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
            'C:/Emod/KinGUII/WorkingDirectory/Voles/Cereals/Cereals6118Fixed IRLS SFO  PTZ.kgo')

   kingraph(Fit,
            'C:/Emod/KinGUII/WorkingDirectory/Voles/Cereals/Cereals6118Fixed IRLS SFO  PTZ.kgg')


#
# End of R-script
#

