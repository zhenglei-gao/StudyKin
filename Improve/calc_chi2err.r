## Calculate chi2 error levels according to FOCUS (2006)
## 0 values at sample time 0 should not be used.
calc_chi2err <- function(observed)
{
  observed1 <- observed
  observed1 <- observed[!(observed$time==0 & observed$value==0),]
  means <- aggregate(value ~ time + name, data = observed1, mean, na.rm=TRUE)##using the mean of repeated measurements.
  ##browser()
  ##errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
  errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
  ## !!!here is the problem!!! observed has two values, thus not be able to really merge!!!!
  ## errdata <- merge(errdata, observed,by = c("time", "name"))
  ## names(errdata)[5] <- 'value_obs'
  errobserved <- merge(observed, predicted_long, by = c("time", "name"), suffixes = c("_obs", "_pred"))
  errdata <- errdata[order(errdata$time, errdata$name), ]
  errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim),errobserved)
  errmin <- data.frame(err.min = errmin.overall$err.min,
                       n.optim = errmin.overall$n.optim, df = errmin.overall$df,
                       err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,
                       EF=errmin.overall$EF,R2=errmin.overall$R2)
  rownames(errmin) <- "All data"
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    errobserved.var <- subset(errobserved, name == obs_var)
    if(outpartri=='default'){
      ##n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f", obs_var,'to', sep="_"), names(parms.optim)))
      n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f",'.*','to',obs_var,sep="_"), names(parms.optim)))
    }
    if(outpartri=='water-sediment'){
      n.k.optim <- length(grep(paste("k_", obs_var, '_',sep=""), names(parms.optim)))
    }
    n.initials.optim <- as.numeric((paste('M0_',obs_var, sep="")) %in% (names(state.ini.optim)))#n.initials.optim <- length(grep(paste('M0_',obs_var, sep=""), names(state.ini.optim)))
    n.optim <- n.k.optim + n.initials.optim
    ## ## added
    k1name <- paste("k1", obs_var,  sep="_")
    k2name <- paste("k2", obs_var,  sep="_")
    gname <- paste("g", obs_var,  sep="_")
    tbname <- paste("tb", obs_var,  sep="_")
    alphaname <- paste("alpha", obs_var,  sep="_")
    betaname <- paste("beta", obs_var,  sep="_")
    ## #
    ## if ("alpha" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("beta" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("k1" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("k2" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("g" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("tb" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## #
    if (alphaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (betaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k1name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k2name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (gname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (tbname %in% names(parms.optim)) n.optim <- n.optim + 1
    
    ##                             #errmin.tmp <- mkinerrmin(errdata.var, n.optim)
    errmin.tmp <- chi2err(errdata.var, n.optim,errobserved.var)
    errmin[obs_var, c("err.min", "n.optim", "df",'err.sig','RMSE','EF','R2')] <- errmin.tmp
  }
  fit$errmin <- errmin
  
  
}