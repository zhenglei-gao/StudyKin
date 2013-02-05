
testfit.gui <- function(mkinmodini,
  eigen = FALSE,
  plot = FALSE, plottitle='',quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  ctr=kingui.control(),...)
{
    #### Control parameters ####
    method <- ctr$method
    odesolver <- ctr$odesolver
    atol <- ctr$atol
    rtol <- ctr$rtol
    control <- ctr$control
    ## mkinmodini is an object by mkinmod.gui
    parms.ini <- mkinmodini$parms.ini
    state.ini <- mkinmodini$state.ini
    lower <- mkinmodini$lower
    upper <- mkinmodini$upper
    fixed_parms <- mkinmodini$fixed_parms
    fixed_initials <- mkinmodini$fixed_initials
    mod_vars <- names(mkinmodini$diffs)
    observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
  # Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmodini$map))
  # Get names of observed variables
  obs_vars = unique(as.character(observed$name))

  # Name the parameters if they are not named yet
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms

  # Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Parameters to be optimised
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]
  ### Temprary For comparison with mkinfit

    ## parms.optim <- c(9.163040e-02,    3.756059e-01,    8.669980e+00,    4.090553e-18 ,1.766376e-02,    1.164488e-02,   3.314102e-01,    6.298495e-01,    1.484640e-01,1.184215e-01,    1.729477e-05,    9.972716e-01,    2.134810e-02,    1.976447e-02 )

    ## names(parms.optim) <- optim_parms
    #print(parms.optim)
     ### ### ### ### ### ###
  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  state.ini.optim.boxnames <- names(state.ini.optim)
    state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
      names(state.ini.optim) <- paste('M0',names(state.ini.optim),  sep="_")
  }
    if(length(state.ini.fixed) > 0) {
      names(state.ini.fixed) <- paste('M0',names(state.ini.fixed), sep="_")
  }
    ### temporary
    ## state.ini.optim <- 1.068602e+02
    ## if(length(state.ini.optim) > 0) {
    ##     names(state.ini.optim) <- paste(state.ini.optim.boxnames, "0", sep="_")
    ## }
    ## parms.optim <- c(9.163040e-02,    3.756059e-01,    8.669980e+00,    4.090553e-18 ,1.766376e-02,    1.164488e-02,   3.314102e-01,    6.298495e-01,    1.484640e-01,1.184215e-01,    1.729477e-05,    9.972716e-01,    2.134810e-02,    1.976447e-02 )
    ## names(parms.optim) <- optim_parms
    ## state.ini.optim <- 1.068602e+02
    ## if(length(state.ini.optim) > 0) {
    ##     names(state.ini.optim) <- paste("M0",state.ini.optim.boxnames, sep="_")
    ## }

    ## ##### ### #####
  # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }

  # Create a function calculating the differentials specified by the model
  # if necessary
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

   # lower <- mkinmodini$lower
   # upper <- mkinmodini$upper
  cost.old <- 1e100
  calls <- 0
  out_predicted <- NA

  # Define the model cost function
  cost <- function(P)
  {
      #names(P) <- pnames
    assign("calls", calls+1, inherits=TRUE)
    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames,state.ini.fixed.boxnames)
    } else {
        odeini <- state.ini.fixed
        names(odeini) <- c( state.ini.fixed.boxnames)
    }

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    outtimes = unique(observed$time)
    evalparse <- function(string)
    {
      eval(parse(text=string), as.list(c(odeparms, odeini)))
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
      ## #sink("NUL")
      ## if(diagnostics(out)$istate[1]!=2)
      ##     {
      ##         #browser()
      ##          out <- ode(
      ##                     y = odeini,
      ##                     times = outtimes,
      ##                     func = mkindiff,
      ##                     parms = odeparms,
      ##                     atol = atol,
      ##                     rtol = rtol,
      ##                     method='ode45'
      ##                     )

      ##     }
      ## #sink()
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
        #browser()
        warning('Integration not completed')
        out_transformed <- apply(out_transformed,2,function(x) {if(sum(is.nan(x))>nrow(out_transformed)-2) x <- rep(Inf,nrow(out_transformed)) else x <- x})
    }
    if(nrow(out_transformed)<length(outtimes))
    {
        tmp <- matrix(0,length(outtimes),ncol(out_transformed)-1)
        tmpnames <-names(out_transformed)
        out_transformed <- data.frame(time=outtimes,tmp)
        names(out_transformed) <- tmpnames
    }

    mC <- modCost(out_transformed, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (cost.old-mC$model > 1) {
      if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

      # Plot the data and current model output if requested
      if(plot) {
        outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)
        if (solution == "analytical") {
          o_plot <- switch(parent.type,
            SFO = SFO.solution(outtimes_plot,
                evalparse(parent.name),
                evalparse(paste("k", parent.name,  sep="_"))),
            FOMC = FOMC.solution(outtimes_plot,
                evalparse(parent.name),
                evalparse(paste("alpha", parent.name,  sep="_")),evalparse(paste("beta", parent.name,  sep="_"))),
            DFOP = DFOP.solution(outtimes_plot,
                evalparse(parent.name),
                evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("g", parent.name,  sep="_"))),
            HS = HS.solution(outtimes_plot,
                evalparse(parent.name),
                evalparse(paste("k1", parent.name,  sep="_")),evalparse(paste("k2", parent.name,  sep="_")),evalparse(paste("tb", parent.name,  sep="_"))),
            SFORB = SFORB.solution(outtimes_plot,
                evalparse(parent.name),
                evalparse(paste("k", parent.name, "bound", sep="_")),
                evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
                evalparse(paste("k", parent.name, "sink", sep="_")))
          )
          out_plot <- cbind(outtimes_plot, o_plot)
          dimnames(out_plot) <- list(outtimes_plot, c("time", sub("_free", "", parent.name)))
      }
        if(solution == "eigen") {
          o_plot <- matrix(mapply(f.out, outtimes_plot),
            nrow = length(mod_vars), ncol = length(outtimes_plot))
          dimnames(o_plot) <- list(mod_vars, outtimes_plot)
          out_plot <- cbind(time = outtimes_plot, t(o_plot))
        }
        if (solution == "deSolve") {
          out_plot <- ode(
            y = odeini,
            times = outtimes_plot,
            func = mkindiff,
            parms = odeparms)
        }
        out_transformed_plot <- data.frame(time = out_plot[,"time"])
        for (var in names(mkinmodini$map)) {
          if((length(mkinmodini$map[[var]]) == 1) || solution == "analytical") {
            out_transformed_plot[var] <- out_plot[, var]
          } else {
            out_transformed_plot[var] <- rowSums(out_plot[, mkinmodini$map[[var]]])
          }
        }

        plot(0, type="n",
          xlim = range(observed$time), ylim = range(observed$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed",main=plottitle)
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)),
            pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_transformed_plot$time, out_transformed_plot[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars,
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
    }

      assign("cost.old", mC$model, inherits=TRUE)
  }
    return(mC)
}
    if(plot) x11()
    ## if(method=='trust') {
    ##     outtimes <-unique(observed$time)
    ##     fit <- modFit1(cost1, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method=method,control=control, state.ini.optim,state.ini.fixed,state.ini.optim.boxnames, state.ini.fixed.boxnames,parms.fixed,outtimes,mkindiff,mkinmodini,observed,err,weight,scaleVar) }else  fit <- modFit1(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method=method,control=control,...)
    #browser()
    if(method=='BFGS')
    {
        parms.optim['k_Parent'] <- 0.1
        fit <- modFit1(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method='BFGS',control=control,...)
    }
##     if(method=='DEoptim')
##     {
##     pnames=names(c(state.ini.optim, parms.optim))
##     fn <- function(P,...){
##         names(P) <- pnames
##         FF<<-cost(P,...)
##         return(FF$model)}
##     upper=c(300,1000,1000,1)
##     fit <- DEoptim(fn, lower, upper)
## }
##     if(method=='Rvmmin')
##     {
##         pnames=names(c(state.ini.optim, parms.optim))
##         fn <- function(P,...){
##             names(P) <- pnames
##             FF<<-cost(P,...)
##             return(FF$model)}
##         upper=c(300,1000,1000,1)
##         fit <- Rvmmin(c(state.ini.optim, parms.optim),fn,lower,upper ,control=list(usenumDeriv=TRUE))
## }
   if(method=='solnp')
{
    pnames=names(c(state.ini.optim, parms.optim))
    fn <- function(P){
        names(P) <- pnames
        FF<<-cost(P)
        return(FF$model)}
    a <- try(fit <- solnp(c(state.ini.optim, parms.optim),fun=fn,LB=lower,UB=upper,control=control),silent=TRUE)
    if(class(a) == "try-error")
    {
        print('solnp fails, try hee other algorithm by users choice, might take longer time. Do something else!')
        warning('solnp fails, switch to  PORT or other algorithm by users choice')
        fit <- modFit1(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method='Marq',control=list(maxIter=200,ftol=1e-9))
         a <- try(fit <- solnp(fit$par,fun=fn,LB=lower,UB=upper,control=control),silent=TRUE)
    }else{
        ### other list need to be attached to fit to give comparable results as in modFit.
        fit$ssr <- fit$values[length(fit$values)]
        fit$residuals <-FF$residual$res
        # mean square per varaible
        if (class(FF) == "modCost") {
            names(fit$residuals)  <- FF$residuals$name
            fit$var_ms            <- FF$var$SSR/FF$var$N
            fit$var_ms_unscaled   <- FF$var$SSR.unscaled/FF$var$N
            fit$var_ms_unweighted <- FF$var$SSR.unweighted/FF$var$N

            names(fit$var_ms_unweighted) <- names(fit$var_ms_unscaled) <-
                names(fit$var_ms) <- FF$var$name
        } else fit$var_ms <- fit$var_ms_unweighted <- fit$var_ms_unscaled <- NA
        np <- length(c(state.ini.optim, parms.optim))
        fit$rank <- np
        fit$df.residual <- length(fit$residuals) - fit$rank
        ########### Calculating the unscaled covariance ###########
        covar <- try(solve(0.5*fit$hessian), silent = TRUE)   # unscaled covariance
        if(!is.numeric(covar)){
            message <- "Cannot estimate covariance directly from hessian of the optimization"
            warning(message)
            fit$solnp.hessian <- fit$hessian
            jac <- NULL
            if (! is.null(jac))Jac <- jac(res$par)else Jac <- gradient(fn, fit$par, centered = TRUE, ...)
            fit$hessian <- 2 * t(Jac) %*% Jac
            covar <- try(solve(0.5*fit$hessian), silent = TRUE)
             if(!is.numeric(covar)){
                 message <- "Cannot estimate covariance  from hessian calculated by gradient"
                 warning(message)
                 fit$Jac <- Jac
                 covar <- matrix(data = NA, nrow = np, ncol = np)
             }

        }else{
            message <- "ok"
        }
        rownames(covar) <- colnames(covar) <-pnames
        fit$covar <- covar
    }

}

    ## Other choices in nloptr
    ## upper=c(Inf,Inf,Inf,1)
    ## fit <- nlm(fn,c(state.ini.optim, parms.optim) )
    ## parms.optim['k_Parent'] <- 0.1
    ## system.time(fit <- nloptr(c(state.ini.optim, parms.optim),eval_f=fn,lb=lower,ub=upper,opts=list(algorithm='NLOPT_GN_DIRECT',maxeval=1000*100,local_opts=list(algorithm='NLOPT_LN_NELDERMEAD'))))
    ##  system.time(fit <- nloptr(c(state.ini.optim, parms.optim),eval_f=fn,lb=lower,ub=upper,opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=1000*100)))


 # We need to return some more data for summary and plotting
  fit$solution <- solution
  if (solution == "eigen") {
    fit$coefmat <- mkinmodini$coefmat
  }
  if (solution == "deSolve") {
    fit$mkindiff <- mkindiff
  }

  # We also need various other information for summary and plotting
  fit$map <- mkinmodini$map
  fit$diffs <- mkinmodini$diffs
  fit$observed <- mkinmodini$residue
   if(method=='trust'){
       P <- fit$par
        if (length(state.ini.optim) > 0) {
            odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
            names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
        }
        else odeini <- state.ini.fixed
        odeparms <- c(P[(length(state.ini.optim) + 1):length(P)],
            parms.fixed)
        #outtimes = unique(observed$time)
        out <- ode(y = odeini, times = outtimes, func = mkindiff,
            parms = odeparms)
        out_transformed <- data.frame(time = out[, "time"])
        for (var in names(mkinmodini$map)) {
            if (length(mkinmodini$map[[var]]) == 1) {
                out_transformed[var] <- out[, var]
            }
            else {
                out_transformed[var] <- rowSums(out[, mkinmodini$map[[var]]])
            }
        }
        assign("out_predicted", out_transformed, inherits = TRUE)
   }
    predicted_long <- mkin_wide_to_long(out_predicted, time = "time")
    fit$predicted <- out_predicted

  # Collect initial parameter values in two dataframes
  fit$start <- data.frame(initial = c(state.ini.optim, parms.optim))
  fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
  fit$start$lower <- lower
  fit$start$upper <- upper

  fit$fixed <- data.frame(
    value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))

  # Calculate chi2 error levels according to FOCUS (2006)
  means <- aggregate(value ~ time + name, data = observed, mean, na.rm=TRUE)##using the mean of repeated measurements.
    #browser()
  #errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
    errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
    errdata <- merge(errdata, observed,by = c("time", "name"))
    names(errdata)[5] <- 'value_obs'
  errdata <- errdata[order(errdata$time, errdata$name), ]
  errmin.overall <- chi2err(errdata, length(parms.optim) + length(state.ini.optim))
  errmin <- data.frame(err.min = errmin.overall$err.min,
    n.optim = errmin.overall$n.optim, df = errmin.overall$df,err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,EF=errmin.overall$EF,R2=errmin.overall$R2)
  rownames(errmin) <- "All data"
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    n.k.optim <- length(grep(paste("k", obs_var, sep="_"), names(parms.optim)))+length(grep(paste("f", obs_var, sep="__"), names(parms.optim)))
    n.initials.optim <- length(grep(paste(obs_var, ".*", "_0", sep=""), names(state.ini.optim)))
    n.optim <- n.k.optim + n.initials.optim
    #### added
    k1name <- paste("k1", obs_var,  sep="_")
    k2name <- paste("k2", obs_var,  sep="_")
    gname <- paste("g", obs_var,  sep="_")
    tbname <- paste("tb", obs_var,  sep="_")
    alphaname <- paste("alpha", obs_var,  sep="_")
    betaname <- paste("beta", obs_var,  sep="_")
    ###
    ## if ("alpha" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("beta" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("k1" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("k2" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("g" %in% names(parms.optim)) n.optim <- n.optim + 1
    ## if ("tb" %in% names(parms.optim)) n.optim <- n.optim + 1
    ###
    if (alphaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (betaname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k1name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (k2name %in% names(parms.optim)) n.optim <- n.optim + 1
    if (gname %in% names(parms.optim)) n.optim <- n.optim + 1
    if (tbname %in% names(parms.optim)) n.optim <- n.optim + 1

    #errmin.tmp <- mkinerrmin(errdata.var, n.optim)
    errmin.tmp <- chi2err(errdata.var, n.optim)
    errmin[obs_var, c("err.min", "n.optim", "df",'err.sig','RMSE','EF','R2')] <- errmin.tmp
}
  fit$errmin <- errmin

  # Calculate dissipation times DT50 and DT90 and formation fractions
    parms.all = c(fit$par, parms.fixed)
    fit$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), DT90 = rep(NA, length(obs_vars)),Kinetic=rep(NA,length(obs_vars)),row.names = obs_vars)
    fit$ff <- vector()
    ff_names = names(mkinmodini$ff)
      for (ff_name in ff_names)
      {
        fit$ff[[ff_name]] =
          eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
      }
    #####
  for (obs_var in obs_vars) {
      f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
      fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])
      type = names(mkinmodini$map[[obs_var]])[1]
      k1name <- paste("k1", obs_var,  sep="_")
      k2name <- paste("k2", obs_var,  sep="_")
      gname <- paste("g", obs_var,  sep="_")
      tbname <- paste("tb", obs_var,  sep="_")
      alphaname <- paste("alpha", obs_var,  sep="_")
      betaname <- paste("beta", obs_var,  sep="_")
    if (type == "SFO") {
      #k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
      #k_tot = sum(parms.all[k_names])
        k_name <- grep(paste("k", obs_var,sep="_"), names(parms.all), value=TRUE)
        k_tot <- parms.all[k_name]
      DT50 = log(2)/k_tot
      DT90 = log(10)/k_tot
      ## for (k_name in k_names)
      ## {
      ##   fit$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
      ## }
    }
    if (type == "FOMC") {
      ## alpha = parms.all["alpha"]
      ## beta = parms.all["beta"]
      alpha = parms.all[alphaname]
      beta = parms.all[betaname]
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      ## ff_names = names(mkinmodini$ff)
      ## for (ff_name in ff_names)
      ## {
      ##   fit$ff[[paste(obs_var, ff_name, sep="_")]] =
      ##     eval(parse(text = mkinmodini$ff[ff_name]), as.list(parms.all))
      ## }
      ## fit$ff[[paste(obs_var, "sink", sep="_")]] = 1 - sum(fit$ff)
    }
    if (type == "DFOP") {
      ## k1 = parms.all["k1"]
      ## k2 = parms.all["k2"]
      ## g = parms.all["g"]
      k1 = parms.all[k1name]
      k2 = parms.all[k2name]
      g = parms.all[gname]
      f <- function(t, x) {
        ((g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)) - (1 - x/100))^2
      }
    }
    if (type == "HS") {
      ## k1 = parms.all["k1"]
      ## k2 = parms.all["k2"]
      ## tb = parms.all["tb"]
      k1 = parms.all[k1name]
      k2 = parms.all[k2name]
      tb = parms.all[tbname]
      f <- function(t, x) {
	fraction = ifelse(t <= tb, exp(-k1 * t), exp(-k1 * tb) * exp(-k2 * (t - tb)))
	(fraction - (1 - x/100))^2
      }
      ##DT50=1
      ##DT90=2
    }
    if (type %in% c("DFOP", "HS")) {
      DTmax <- 1000
      DT50.o <- optimize(f, c(0.001, DTmax), x=50)$minimum
      DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
      DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
      DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)

    }
    if (type == "SFORB") {
      # FOCUS kinetics (2006), p. 60 f
      k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
      k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
      k_1output = sum(parms.all[k_out_names])
      k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
      k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]

      sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
      b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
      b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

      SFORB_fraction = function(t) {
        ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
      }
      f_50 <- function(t) (SFORB_fraction(t) - 0.5)^2
      max_DT <- 1000
      DT50.o <- optimize(f_50, c(0.01, max_DT))$minimum
      if (abs(DT50.o - max_DT) < 0.01) DT50 = NA else DT50 = DT50.o
      f_90 <- function(t) (SFORB_fraction(t) - 0.1)^2
      DT90.o <- optimize(f_90, c(0.01, 1000))$minimum
      if (abs(DT90.o - max_DT) < 0.01) DT90 = NA else DT90 = DT90.o
      for (k_out_name in k_out_names)
      {
        fit$ff[[sub("k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
      }
  }
    fit$distimes[obs_var, ] = c(DT50, DT90,type)
}
#browser()
  # Collect observed, predicted and residuals
  data <- merge(observed, predicted_long, by = c("time", "name"))
  names(data) <- c("time", "variable", "observed", "predicted")
  data$residual <- data$observed - data$predicted
  data$variable <- ordered(data$variable, levels = obs_vars)
  fit$data <- data[order(data$variable, data$time), ]
  fit$atol <- atol

class(fit) <- c('kingui',"mkinfit", "modFit")
return(fit)
}

## NLOPT_GN_DIRECT
## NLOPT_GN_DIRECT_L
## NLOPT_GN_DIRECT_L_RAND
## NLOPT_GN_DIRECT_NOSCAL
## NLOPT_GN_DIRECT_L_NOSCAL
## NLOPT_GN_DIRECT_L_RAND_NOSCAL
## NLOPT_GN_ORIG_DIRECT
## NLOPT_GN_ORIG_DIRECT_L
## NLOPT_GD_STOGO
## NLOPT_GD_STOGO_RAND
## NLOPT_LD_SLSQP
## NLOPT_LD_LBFGS_NOCEDAL
## NLOPT_LD_LBFGS
## NLOPT_LN_PRAXIS
## NLOPT_LD_VAR1
## NLOPT_LD_VAR2
## NLOPT_LD_TNEWTON
## NLOPT_LD_TNEWTON_RESTART
## NLOPT_LD_TNEWTON_PRECOND
## NLOPT_LD_TNEWTON_PRECOND_RESTART
## NLOPT_GN_CRS2_LM
## NLOPT_GN_MLSL
## NLOPT_GD_MLSL
## NLOPT_GN_MLSL_LDS
## NLOPT_GD_MLSL_LDS
## NLOPT_LD_MMA
## NLOPT_LN_COBYLA
## NLOPT_LN_NEWUOA
## NLOPT_LN_NEWUOA_BOUND
## NLOPT_LN_NELDERMEAD
## NLOPT_LN_SBPLX
## NLOPT_LN_AUGLAG
## NLOPT_LD_AUGLAG
## NLOPT_LN_AUGLAG_EQ
## NLOPT_LD_AUGLAG_EQ
## NLOPT_LN_BOBYQA
## NLOPT_GN_ISRES
