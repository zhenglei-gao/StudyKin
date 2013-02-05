modplot <- function(mkinmodini,eigen=FALSE,name=NULL, xlab = "Time", ylab = "Observed", legend = TRUE,pdf=TRUE,  err = NULL, weight = "none",scaleVar = FALSE,ctr=kingui.control(),...)
{
    ## example usage:
    ## kinplot(mod,name=c('MESO','M851','M459','M460','M095','M944'))
#### Control parameters ####
    xlim <-c(1,1.1)*range(mkinmodini$residue$time)
    ylim = c(1,1.1)*range(mkinmodini$residue[,2:ncol(mkinmodini$residue)], na.rm = TRUE)
    outtimes1 <- seq(xlim[1], xlim[2], length.out=100)
    method <- ctr$method
    odesolver <- ctr$odesolver
    atol <- ctr$atol
    rtol <- ctr$rtol
    control <- ctr$control
### This is a modification based on the "IRLSkinfit0" function.
### version
    parms.ini <- mkinmodini$parms.ini
    state.ini <- mkinmodini$state.ini
    lower <- mkinmodini$lower
    upper <- mkinmodini$upper
    fixed_parms <- mkinmodini$fixed_parms
    fixed_initials <- mkinmodini$fixed_initials

    mod_vars <- names(mkinmodini$diffs)
    observed <- mkin_wide_to_long(mkinmodini$residue,time='time')
    observed <- subset(observed, name %in% names(mkinmodini$map))
    NAind <-which(is.na(observed$value))
    ERR <- rep(1,nrow(observed))
    observed <- cbind(observed,err=ERR)
    obs_vars = unique(as.character(observed$name))
    if (is.null(names(parms.ini)))   names(parms.ini) <- mkinmodini$parms

    if (is.null(names(state.ini)))  names(state.ini) <- mod_vars
    parms.fixed <- parms.ini[fixed_parms]
    optim_parms <- setdiff(names(parms.ini), fixed_parms)
    parms.optim <- parms.ini[optim_parms]
    state.ini.fixed <- state.ini[fixed_initials]
    state.ini.fixed.boxnames <- names(state.ini.fixed)
    optim_initials <- setdiff(names(state.ini), fixed_initials)
    state.ini.optim <- state.ini[optim_initials]
    state.ini.optim.boxnames <- names(state.ini.optim)
    if (length(state.ini.optim) > 0) {
        names(state.ini.optim) <- paste('M0',names(state.ini.optim),sep = "_")
    }
     if(length(state.ini.fixed) > 0) {
      names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }

## ### temporary
##      parms.optim <- c(9.163040e-02,    3.756059e-01,    8.669980e+00,    4.090553e-18 ,1.766376e-02,    1.164488e-02,   3.314102e-01,    6.298495e-01,    1.484640e-01,1.184215e-01,    1.729477e-05,    9.972716e-01,    2.134810e-02,    1.976447e-02 )
##     names(parms.optim) <- optim_parms
##     state.ini.optim <- 1.068602e+02
##     if(length(state.ini.optim) > 0) {
##         names(state.ini.optim) <- paste(state.ini.optim.boxnames, "0", sep="_")
##     }

    # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
   #browser()
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }

 if(solution == "deSolve") {
    mkindiff <- function(t, state, parms) {
      time <- t
      diffs <- vector()
      for (box in mod_vars)
      {
        diffname <- paste("d", box, sep="_")
        diffs[diffname] <- with(as.list(c(time,state, parms)),
          eval(parse(text=mkinmodini$diffs[[box]])))
      }
      return(list(c(diffs)))
    }
  }
    if(length(state.ini.optim) > 0) {
      odeini <- c(state.ini.optim,state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
        odeini <- state.ini.fixed
        names(odeini) <- c( state.ini.fixed.boxnames)
    }

    odeparms <- c(parms.optim, parms.fixed)

    outtimes = unique(observed$time)
    evalparse <- function(string)
    {
      eval(parse(text=string), as.list(c(odeparms, odeini)))
    }

    # Solve the system
    if (solution == "analytical") {
      parent.type = names(mkinmodini$map[[1]])[1]
      parent.name = names(mkinmodini$diffs)[[1]]
       o <- switch(parent.type,
        SFO = SFO.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k", parent.name,  sep="_"))),
        FOMC = FOMC.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("alpha", parent.name,  sep="_")), evalparse(paste("beta", parent.name,  sep="_"))),
        DFOP = DFOP.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k1", parent.name,  sep="_")), evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
        HS = HS.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
        SFORB = SFORB.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k", parent.name, "bound", sep="_")),
            evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
            evalparse(paste("k", parent.name, "sink", sep="_")))
      )


      out <- cbind(outtimes, o)
      dimnames(out) <- list(outtimes, c("time", sub("_free", "", parent.name)))
    }
    if (solution == "eigen") {
      coefmat.num <- matrix(sapply(as.vector(mkinmodini$coefmat), evalparse),
        nrow = length(mod_vars))
      e <- eigen(coefmat.num)
      zz.ev <- e$values
      if(min(zz.ev)[1]<0){

        warning("\'coefmat is not positive definite!\n")
        solution <- 'deSolve' ## switch to deSolve methods
        if(solution == "deSolve") {
            mkindiff <- function(t, state, parms) {
                time <- t
                diffs <- vector()
                for (box in mod_vars)
                {
                    diffname <- paste("d", box, sep="_")
                    diffs[diffname] <- with(as.list(c(time,state, parms)),
                                            eval(parse(text=mkinmodini$diffs[[box]])))
                }
                return(list(c(diffs)))
            }
        }
    }else{
        cc <- solve(e$vectors, odeini)
        f.out <- function(t) {
            e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% cc
        }
        o <- matrix(mapply(f.out, outtimes),
                    nrow = length(mod_vars), ncol = length(outtimes))
        dimnames(o) <- list(mod_vars, outtimes)
        out <- cbind(time = outtimes, t(o))
    }
  }
    if (solution == "deSolve")
    {
      out <- ode(
        y = odeini,
        times = outtimes,
        func = mkindiff,
        parms = odeparms,
        atol = atol,
        rtol = rtol,
        method=odesolver
      )
    }

    # Output transformation for models with unobserved compartments like SFORB
    out_transformed <- data.frame(time = out[,"time"])
    for (var in names(mkinmodini$map)) {
      if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
        out_transformed[var] <- out[, var]
      } else {
        out_transformed[var] <- rowSums(out[, mkinmodini$map[[var]]])
      }
    }
    assign("out_predicted", out_transformed, inherits=TRUE)
    if(sum(apply(out_transformed,2,function(x) sum(is.nan(x))>nrow(out_transformed)-2))>0)
    {
        out_transformed <- apply(out_transformed,2,function(x) {if(sum(is.nan(x))>nrow(out_transformed)-2) x <- rep(Inf,nrow(out_transformed)) else x <- x})
    }
    if(nrow(out_transformed)<length(outtimes))
    {
        tmp <- matrix(0,length(outtimes),ncol(out_transformed)-1)
        tmpnames <-names(out_transformed)
        out_transformed <- data.frame(time=outtimes,tmp)
        names(out_transformed) <- tmpnames
    }
    #browser()
    mC <- modCost(out_transformed, observed, y = "value",
      err = 'err', weight = weight, scaleVar = scaleVar)

    out_predicted <- out_transformed
    predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
    observed <-  mkinmodini$residue
    data <- merge(mkin_wide_to_long(observed,'time'), predicted_long, by = c("time", "name"))
    names(data) <- c("time", "variable", "observed","predicted")
    data$variable <- ordered(data$variable, levels = obs_vars)
    data <- data[order(data$variable, data$time), ]
#########################
outtimes = outtimes1

    # Solve the system
    if (solution == "analytical") {
      parent.type = names(mkinmodini$map[[1]])[1]
      parent.name = names(mkinmodini$diffs)[[1]]
       o <- switch(parent.type,
        SFO = SFO.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k", parent.name,  sep="_"))),
        FOMC = FOMC.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("alpha", parent.name,  sep="_")), evalparse(paste("beta", parent.name,  sep="_"))),
        DFOP = DFOP.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k1", parent.name,  sep="_")), evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
        HS = HS.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
        SFORB = SFORB.solution(outtimes,
            evalparse(parent.name),
            evalparse(paste("k", parent.name, "bound", sep="_")),
            evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
            evalparse(paste("k", parent.name, "sink", sep="_")))
      )


      out <- cbind(outtimes, o)
      dimnames(out) <- list(outtimes, c("time", sub("_free", "", parent.name)))
    }
    if (solution == "eigen") {
      coefmat.num <- matrix(sapply(as.vector(mkinmodini$coefmat), evalparse),
        nrow = length(mod_vars))
      e <- eigen(coefmat.num)
      zz.ev <- e$values
      if(min(zz.ev)[1]<0){

        warning("\'coefmat is not positive definite!\n")
        solution <- 'deSolve' ## switch to deSolve methods
        if(solution == "deSolve") {
            mkindiff <- function(t, state, parms) {
                time <- t
                diffs <- vector()
                for (box in mod_vars)
                {
                    diffname <- paste("d", box, sep="_")
                    diffs[diffname] <- with(as.list(c(time,state, parms)),
                                            eval(parse(text=mkinmodini$diffs[[box]])))
                }
                return(list(c(diffs)))
            }
        }
    }else{
        cc <- solve(e$vectors, odeini)
        f.out <- function(t) {
            e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% cc
        }
        o <- matrix(mapply(f.out, outtimes),
                    nrow = length(mod_vars), ncol = length(outtimes))
        dimnames(o) <- list(mod_vars, outtimes)
        out <- cbind(time = outtimes, t(o))
    }
  }
    if (solution == "deSolve")
    {
      out <- ode(
        y = odeini,
        times = outtimes,
        func = mkindiff,
        parms = odeparms,
        atol = atol,
        rtol = rtol,
        method=odesolver
      )
    }

    # Output transformation for models with unobserved compartments like SFORB
    out_transformed <- data.frame(time = out[,"time"])
    for (var in names(mkinmodini$map)) {
      if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
        out_transformed[var] <- out[, var]
      } else {
        out_transformed[var] <- rowSums(out[, mkinmodini$map[[var]]])
      }
    }
    assign("out_predicted", out_transformed, inherits=TRUE)
    if(sum(apply(out_transformed,2,function(x) sum(is.nan(x))>nrow(out_transformed)-2))>0)
    {
        out_transformed <- apply(out_transformed,2,function(x) {if(sum(is.nan(x))>nrow(out_transformed)-2) x <- rep(Inf,nrow(out_transformed)) else x <- x})
    }
    if(nrow(out_transformed)<length(outtimes))
    {
        tmp <- matrix(0,length(outtimes),ncol(out_transformed)-1)
        tmpnames <-names(out_transformed)
        out_transformed <- data.frame(time=outtimes,tmp)
        names(out_transformed) <- tmpnames
    }

########################
if(is.null(name)){
  # Plot the data and model output
  plot(0, type="n",
    xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, ...)
  col_obs <- pch_obs <- 1:length(mkinmodini$map)
  names(col_obs) <- names(pch_obs) <- names(mkinmodini$map)
  for (obs_var in names(mkinmodini$map)) {
    points(subset(data, variable == obs_var, c(time, observed)),
    pch = pch_obs[obs_var], col = col_obs[obs_var])
  }
  matlines(out_transformed$time, out_transformed[-1])
  if (legend == TRUE) {
    legend("topright", inset=c(0.05, 0.05), legend=names(mkinmodini$map),
      col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
  }
}else{
    col_obs <- pch_obs <- 1:length(mkinmodini$map)
    names(col_obs) <- names(pch_obs) <- names(mkinmodini$map)
    for (obs_var in name)
    {
        if(pdf==F) x11()
        plot(subset(data, variable == obs_var, c(time, observed)),
    pch = pch_obs[obs_var],xlim=xlim,ylim=c(1,1.1)*range(subset(data, variable == obs_var, c( observed,predicted)),na.rm=TRUE), col = col_obs[obs_var],main=obs_var)
        #browser()
        lines(out_transformed$time,as.vector(out_transformed[obs_var][[1]]))
    }
}

return(mC)
}
