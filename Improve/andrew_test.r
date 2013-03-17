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




Test <- mkinmod.full(
  Parent = list(
    time = c(    0,     0,     1,     1,     3,     3,     7,     7,    16,    16,    21,    21,    28,    28),
    residue = c(1.676, 2.095, 1.362, 2.444, 0.628, 0.506, 0.310, 0.235, 0.075, 0.088, 0.030, 0.033, 0.023, 0.015),
    weight = c(    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1),
    to = c("Metab"),
    FF = list(ini   = c(0.50),
              fixed = c(0),
              lower = c(0.0),
              upper = c(1.0)),
    sink  = TRUE,
    type = "SFO",
    k = list(ini   = 0.3000,
             fixed = 0,
             lower = 0.0,
             upper = Inf),
    M0 = list(ini   = 1.676,
              fixed = 0,
              lower = 0.0,
              upper = Inf)),
  Metab = list(
    time = c(    0,     0,     1,     1,     3,     3,     7,     7,    16,    16,    21,    21,    28,    28),
    residue = c(0.192, 0.219, 0.461, 0.444, 0.375, 0.294, 0.170, 0.175, 0.034, 0.051, 0.020, 0.026, 0.014, 0.011),
    weight = c(    1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1),
    sink  = TRUE,
    type = "SFO",
    k = list(ini   = 0.9000,
             fixed = 0,
             lower = 0.0,
             upper = Inf),
    M0 = list(ini   = 0,
              fixed = 0,
              lower = 0.0,
              upper = Inf)))

#
# Fit and optimizer
#

Fit    <- IRLSkinfit.full(
  Test,
  plot      = TRUE,
  quiet     = TRUE,
  ctr       = kingui.control(
    method = 'L-BFGS-B',
    submethod = 'Port',
    maxIter = 100,
    tolerance = 1E-06,
    odesolver = 'lsoda'),
  irls.control = list(
    maxIter = 10,
    tolerance = 0.001))

new1 <- compare_kin_mod(Test,Fit=Fit)$new_mkinmod
a <- compare_kin_mod(new1,newparms=c(2.0633,0.2105,0.3033,0.8945,1))
a$comparison
