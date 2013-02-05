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


sourceDir('D:/Projects/Kinetics/KINGUI/package/Sourcefunctions/')
version <- readversion('D:/Projects/Kinetics/KINGUI/package/Sourcefunctions/')


############
ws <- read.table('ws.txt')
wsmodel <- mkinmod.full(P_water=list(type='SFO',to=c('P_sediment','M_water')),
                        P_sediment=list(type='SFO',to=c('P_water','M_sediment')),
                        M_water=list(type='SFO',to=c('M_sediment')),
                        M_sediment=list(type='SFO',to=c('M_water')),
                        ##inpartri='water-sediment',outpartri='water-sediment',
                         inpartri='default',outpartri='water-sediment',
                        data=ws
                        )
 Fit    <- IRLSkinfit.full(
                           wsmodel,
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

list2ascii(summary(Fit, cov = TRUE,version=version),'water_sediment_example.kgo')
