
#require(Rsolnp)
#require(optimx)
#require(DEoptim)


##' Fit a kinetic model using the IRLS algorithm.
##'
##' This function does kinetic evaluations using the IRLS algorithm.
##' This is deprecated! Use \code{\link{IRLSkinfit.full}} instead.
##' @aliases IRLSkinfit.gui summary.kingui
##' @param mkinmodini A list of class
##' \code{\link{mkinmod.gui}}, containing the kinetic model to be fitted to the
##' data, and the initial parameter values, the observed data.
##' @param eigen If TRUE, the solution of the
##' system of differential equations should be based on the spectral
##' decomposition of the coefficient matrix in cases that this is possible.
##' @param plot If TRUE,the observed values
##' and the numerical solutions should be plotted at each stage of the
##' optimisation.
##' @param plottitle The title of the
##' plot for visualizing the optimization process.
##' @param quiet If TRUE, suppress printing
##' out the current model cost after each(>1) improvement.
##' @param err See argumetns of
##' \code{\link{mkinfit.gui}}
##' @param weight See argumetns of
##' \code{\link{mkinfit.gui}}
##' @param scaleVar See argumetns of
##' \code{\link{mkinfit.gui}}
##' @param ctr a list of control values for the
##' estimation algorithm to replace the default values including maximum
##' iterations and absolute error tolerance.  Defaults to the output of
##' \code{\link{kingui.control}}.
##' @param irls.control  A list of
##' control values for the estimation algorithm to replace the default values
##' including the maximum number of iterations for the outer iteration and the
##' error tolerance level for the error variance estimation updating.
##' @param update %% ~~Describe \code{update} here~~ If not NULL, should be a
##' list of starting values obtained from other optimization methods.
##' @param \dots Further arguments that will
##' be passed to \code{\link{modFit1}}.
##' @return  A list with "kingui", "mkinfit" and "modFit" in the class
##' attribute. A summary can be obtained by \code{\link{summary.kingui}}.
##' @note Deprecated!
##' @author  Zhenglei Gao
##' @seealso
##' \code{\link{IRLSkinfit}}, \code{\link{mkinfit.gui}}
##' @keywords Kinetic-Evaluations
##' @examples
##'
##' SFO_SFO_gui <- mkinmod.gui(Parent = list(type = "SFO", to = "Metab", sink = TRUE,
##'                            k = list(ini = 0.1,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 195,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'                            FF = list(ini = c(.1),
##'                       fixed = c(0),
##'                       lower = c(0),
##'                       upper = c(1)),
##'                       time=c(0.0,2.8,   6.2,  12.0,  29.2,  66.8,  99.8,
##' 127.5, 154.4, 229.9, 272.3, 288.1, 322.9),
##'                     residue = c( 157.3, 206.3, 181.4, 223.0, 163.2,
##' 144.7,  85.0,  76.5,  76.4,  51.5,  45.5,  47.3, 42.7)),
##'                            Metab = list(type = "SFO",
##'                            k = list(ini = 0.1   ,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 0,
##'                       fixed = 1,
##'                       lower = 0,
##'                       upper = Inf),
##'                     residue =c( 0.0,  0.0,  0.0,  1.6,  4.0, 12.3, 13.5,
##' 12.7, 11.4, 11.6, 10.9,  9.5,  7.6))                           )
##' fit <- IRLSkinfit.gui(SFO_SFO_gui,plot=T,ctr=kingui.control(maxIter=100,
##' tolerance=1e-6,odesolver='lsoda'))
##'
IRLSkinfit.gui <- function(mkinmodini,
  eigen = FALSE,
  plot = FALSE, plottitle='',quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  ctr=kingui.control(),irls.control=list(),update=NULL,useHsolnp=FALSE,...)
{
    options(warn=-1)
    #### Control parameters ####
    method <- ctr$method
    odesolver <- ctr$odesolver
    atol <- ctr$atol
    rtol <- ctr$rtol
    control <- ctr$control
    marqctr <- ctr$marqctr
    goMarq <- ctr$goMarq
    submethod <- ctr$submethod
    Hmethod1 <- ctr$Hmethod1
    Hmethod2 <- ctr$Hmethod2
    ## mkinmodini is an object by mkinmod.gui
    parms.ini <- mkinmodini$parms.ini
    state.ini <- mkinmodini$state.ini
    lower <- mkinmodini$lower
    upper <- mkinmodini$upper
    fixed_parms <- mkinmodini$fixed_parms
    fixed_initials <- mkinmodini$fixed_initials
    mod_vars <- names(mkinmodini$diffs)
    observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
    observed$err <-c(as.matrix(mkinmodini$weightmat))
  # Subset dataframe with mapped (modelled) variables
    observed <- subset(observed, name %in% names(mkinmodini$map))
    #ERR <- rep(1,nrow(observed))### start weigting 1
    #observed <- cbind(observed,err=ERR)
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
     ### If updating from previous fit.####
    if(!is.null(update))
    {
        parms.optim <- update$par[optim_parms]
        state.ini.optim <-update$par[names(state.ini.optim)]
        ctr <- update$ctr
        ## Control parameters ####
        method <- ctr$method
        odesolver <- ctr$odesolver
        atol <- ctr$atol
        rtol <- ctr$rtol
        control <- ctr$control
    }
 ######################################################################
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
    ## has to change the odeini order since it is different from the mod_vars order.
    odeini <- odeini[mod_vars]
    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

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
      err = 'err', weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (cost.old-mC$model > ctr$quiet.tol) {
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
###############################
    method <- 'solnp'
   if(method=='solnp')
{
    pnames=names(c(state.ini.optim, parms.optim))
    fn <- function(P){
        names(P) <- pnames
        FF<<-cost(P)
        return(FF$model)}
    a <- try(fit <- solnp(c(state.ini.optim, parms.optim),fun=fn,LB=lower,UB=upper,control=control),silent=TRUE)
    flag <- 1
    if(class(a) == "try-error")
    {
        print('solnp fails, try PORT or other algorithm by users choice, might take longer time. Do something else!')
        warning('solnp fails, switch to  PORT or other algorithm by users choice')
        fit <- modFit1(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, method=submethod,control=kingui.control(method=submethod,tolerance=ctr$control$tol)$control)
        flag <- 0
        ## now using submethod already

    }
        ###########################
        if(length(irls.control)==0) irls.control <- list(maxIter=5,tol=1e-05)
    if(is.null(irls.control$tol)) tol <- 1e-05 else tol <- irls.control$tol
    if(is.null(irls.control$maxIter)) maxIter <- 5 else  maxIter <- irls.control$maxIter



        niter <- 1

        ## insure one IRLS iteration
        diffsigma <- 100
        olderr <- rep(1,length(mod_vars))

        while(diffsigma>tol & niter<=maxIter)
        {
            ### other list need to be attached to fit to give comparable results as in modFit.
            if(flag==1)## fit from solnp
            {
                fit$ssr <- fit$values[length(fit$values)]
                fit$residuals <-FF$residual$res
                ## mean square per varaible
                if (class(FF) == "modCost") {
                    names(fit$residuals)  <- FF$residuals$name
                    fit$var_ms            <- FF$var$SSR/FF$var$N
                    fit$var_ms_unscaled   <- FF$var$SSR.unscaled/FF$var$N
                    fit$var_ms_unweighted <- FF$var$SSR.unweighted/FF$var$N

                    names(fit$var_ms_unweighted) <- names(fit$var_ms_unscaled) <-
                        names(fit$var_ms) <- FF$var$name
                } else fit$var_ms <- fit$var_ms_unweighted <- fit$var_ms_unscaled <- NA
            }
            err1 <- sqrt(fit$var_ms_unweighted)
            ERR <- err1[as.character(observed$name)]
            observed$err <-ERR
            diffsigma <- sum((err1-olderr)^2)
            cat("IRLS iteration at",niter, "; Diff in error variance ", diffsigma,"\n")
            olderr <- err1
            ## #
            if(goMarq==1) {

                print('do a local optmization')
                fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method='Marq',control=marqctr)
                flag <- 0
            }else{
                flag <- 1
                a <- try(fit <- solnp(fit$par,fun=fn,LB=lower,UB=upper,control=control),silent=TRUE)
                if(class(a) == "try-error")
                {
                    flag <- 0
                    print('solnp fails during IRLS iteration, try PORT or other algorithm by users choice, might take longer time. Do something else!')
                    warning('solnp fails during IRLS iteration, switch to  PORT or other algorithm by users choice')

                    fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method=submethod,control=list())
                }
            }
                niter <- niter+1

                ## # If not converged, reweight and fit
        }
    ###########################################
        ### other list need to be attached to fit to give comparable results as in modFit.
        if(flag==1){
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
    }
    ##browser()
        ########### Calculating the unscaled covariance ###########
        if(flag!=1) covar <- try(solve(0.5*fit$hessian), silent = TRUE) else {
            if(useHsolnp==TRUE) covar <- try(solve(0.5*fit$hessian)) else covar <- NULL# unscaled covariance
            fit$solnp.hessian <- fit$hessian
        }
        if(!is.numeric(covar)){
            message <- "Cannot estimate covariance directly from hessian of the optimization"
            warning(message)
            print('Solnp algorithm does not provide a valid hessian. We need to calculate the hessian to estimate the correlation matrix. This may take a while. Please be patient!')
             if(!is.numeric(covar)){
                 fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method=Hmethod1,control=list())
                 covar <- fit$covar
                 if(!is.numeric(covar))
                 {
                       message <- "Cannot estimate covariance  from hessian calculated by gradient or by Hmethod1"
                      warning(message)
                       #print('go to the third level to calculate covar')
                       fit <- modFit1(cost, fit$par, lower = lower, upper = upper, method=Hmethod2,control=list())
                       covar <- fit$covar
                       if(!is.numeric(fit$covar)){
                           covar <- fit$covar
                       }else{
                           covar <- matrix(data = NA, nrow = np, ncol = np)
                           warning('covar not estimable')
                       }
                 }else{
                     ##
             }
             }

        }else{
            message <- "ok"
        }
        rownames(covar) <- colnames(covar) <-pnames
        fit$covar <- covar
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
if(length(state.ini.optim)>0){
    fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=lower[1:length(state.ini.optim)],upper=upper[1:length(state.ini.optim)])}else{
        fit$start0 <- data.frame(initial=state.ini.optim,type=rep("state", length(state.ini.optim)),lower=numeric(0),upper=numeric(0))
    }

        start0 <- mkinmodini$start[mkinmodini$start$fixed==0,]
        fit$start0 <- rbind(fit$start0,data.frame(initial=start0$initial,type=start0$type,lower=start0$lower,upper=start0$upper,row.names =rownames(start0)))
        fit$fixed0 <- data.frame(value = state.ini.fixed,type=rep("state", length(state.ini.fixed)),by=rep("user", length(state.ini.fixed)))
        fixed0 <- mkinmodini$start[mkinmodini$start$fixed==1,]
        if(nrow(fixed0)>0) fit$fixed0 <- rbind(fit$fixed0,data.frame(value=fixed0$initial,type=fixed0$type,by=rep('user',nrow(fixed0)),row.names=rownames(fixed0)))
## keep for the summary part to add the kingui fixed parameters into the optimised parameter section.
  fit$start <- data.frame(initial = c(state.ini.optim, parms.optim))
  fit$start$type = c(rep("state", length(state.ini.optim)), rep("deparm", length(parms.optim)))
  fit$start$lower <- lower
  fit$start$upper <- upper

  fit$fixed <- data.frame(
    value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)), rep("deparm", length(parms.fixed)))
fit$fixed$by <- c(rep("user", length(state.ini.fixed)), mkinmodini$fixed_flag)
#######################################################
  # Calculate chi2 error levels according to FOCUS (2006)
  means <- aggregate(value ~ time + name, data = observed, mean, na.rm=TRUE)##using the mean of repeated measurements.
    #browser()
  #errdata <- merge(means, predicted_long, observed,by = c("time", "name"), suffixes = c("_mean", "_pred",'_obs'))
    errdata <- merge(means, predicted_long, by = c("time", "name"), suffixes = c("_mean", "_pred"))
    errdata <- merge(errdata, observed,by = c("time", "name"))
    names(errdata)[5] <- 'value_obs'
  errdata <- errdata[order(errdata$time, errdata$name), ]
  errmin.overall <- chi2err.bug(errdata, length(parms.optim) + length(state.ini.optim))
  errmin <- data.frame(err.min = errmin.overall$err.min,
    n.optim = errmin.overall$n.optim, df = errmin.overall$df,err.sig = errmin.overall$err.sig,RMSE=errmin.overall$RMSE,EF=errmin.overall$EF,R2=errmin.overall$R2)
  rownames(errmin) <- "All data"
  for (obs_var in obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)
    n.k.optim <- (paste("k", obs_var, sep="_")) %in% (names(parms.optim))+length(grep(paste("f", obs_var,'to', sep="_"), names(parms.optim)))#n.k.optim <- length(grep(paste("k", obs_var, sep="_"), names(parms.optim)))+length(grep(paste("f", obs_var, sep="_"), names(parms.optim)))
    n.initials.optim <- as.numeric((paste('M0_',obs_var, sep="")) %in% (names(state.ini.optim)))#n.initials.optim <- length(grep(paste('M0_',obs_var, sep=""), names(state.ini.optim)))
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
    errmin.tmp <- chi2err.bug(errdata.var, n.optim)
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
##browser()
  for (obs_var in obs_vars) {
      f_tot <- grep(paste(obs_var, "_",sep=''), names(fit$ff), value=TRUE)
      f_exp <- grep(paste(obs_var, "to",obs_var,sep='_'), names(fit$ff), value=TRUE)
      f_exp1 <- grep(paste(obs_var, "to",'sink',sep='_'), names(fit$ff), value=TRUE)
      fit$ff[[paste(obs_var,'to', "sink", sep="_")]] = 1 - sum(fit$ff[f_tot])+sum(fit$ff[f_exp])+sum(fit$ff[f_exp1])
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
        k_name <- paste("k", obs_var,sep="_")#grep(paste("k", obs_var,sep="_"), names(parms.all), value=TRUE)
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
      DTmax <- 1000
      DT50.o <- optimize(f, c(0.0001,DTmax), x=50)$minimum
      DTmax1 <- log(2)/min(k1,k2)
      DT50.o1 <- optimize(f, c(0, DTmax1), x=50)$minimum
      DT50.o <- ifelse(f(DT50.o,50)>f(DT50.o1,50), DT50.o1,DT50.o)
      DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
      DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
      DT90.o1 <- optimize(f, c(0, DTmax1), x=90)$minimum
      DT90.o <- ifelse(f(DT90.o,90)>f(DT90.o1,90), DT90.o1,DT90.o)
      DT90 = ifelse(DTmax - DT90.o < 0.1, NA, DT90.o)
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
      DTmax <- 1000
      hso1 <- optimize(f, c(0.0001,tb), x=50)
      hso2 <- optimize(f, c(tb,DTmax), x=50)
      DT50.o <- ifelse(hso1$objective<=hso2$objective,hso1$minimum,hso2$minimum)
      ## DT50.o <- optimize(f, c(0.0001,tb), x=50)$minimum
      DT50 = ifelse(DTmax - DT50.o < 0.1, NA, DT50.o)
      hso1 <- optimize(f, c(0.0001,tb), x=90)
      hso2 <- optimize(f, c(tb,DTmax), x=90)
      DT90.o <- ifelse(hso1$objective<=hso2$objective,hso1$minimum,hso2$minimum)
      ##DT90.o <- optimize(f, c(0.001, DTmax), x=90)$minimum
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
    fit$distimes[obs_var, ] = c(ifelse(is.na(DT50),NA,formatC(DT50,4,format='f')), ifelse(is.na(DT90),NA,formatC(DT90,4,format='f')),type)
  }
#browser()
  # Collect observed, predicted and residuals
  data <- merge(observed, predicted_long, by = c("time", "name"))
  #names(data) <- c("time", "variable", "observed", "predicted")
  names(data) <- c("time", "variable", "observed","err-std", "predicted")
  data$residual <- data$observed - data$predicted
  data$variable <- ordered(data$variable, levels = obs_vars)
  fit$data <- data[order(data$variable, data$time), ]
  fit$atol <- atol

class(fit) <- c('kingui0',"mkinfit", "modFit")
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


##' S3 method for class 'kingui0'
##'
##'
##' @param object An object of class
##' 'kingui' from the result of NLS or IRLS fit.
##' @param data If TRUE, include in the
##' returned values a data frame containing the observed and predicted values
##' with residuals and estimated standard deviations or weights.
##' @param distimes If TRUE, DT50 and DT90
##' values should be included.
##' @param ff If TRUE, the formation fraction
##' should be calculated from the estimated transformed parameters.
##' @param cov  If TRUE, parameter covariances
##' should be calculated.
##' @param \dots Optional arguments passed to
##' methods like 'print'
##' @return The summary function returns a list with the same
##' components as 'summary.mkinfit', and the additional components from the
##' results of the optimization routine used. %%\item{diffs}{} %%\item{}{}
##' %%\item{}{} %%\item{}{}
##' @note This is deprecated!
##' @author Zhenglei Gao
##' @seealso \code{\link{summary.mkinfit}}, \code{\link{summary.modFit}}
##' @keywords internal
##' @examples
##'
##' \dontrun{
##' SFO_SFO_gui <- mkinmod.gui(Parent = list(type = "SFO", to = "Metab", sink = TRUE,
##'                            k = list(ini = 0.1,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 195,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'                            FF = list(ini = c(.1),
##'                       fixed = c(0),
##'                       lower = c(0),
##'                       upper = c(1)),
##'                       time=c(0.0,2.8,   6.2,  12.0,  29.2,  66.8,  99.8,
##' 127.5, 154.4, 229.9, 272.3, 288.1, 322.9),
##'                     residue = c( 157.3, 206.3, 181.4, 223.0, 163.2,
##' 144.7,  85.0,  76.5,  76.4,  51.5,  45.5,  47.3, 42.7)),
##'                            Metab = list(type = "SFO",
##'                            k = list(ini = 0.1   ,
##'                       fixed = 0,
##'                       lower = 0,
##'                       upper = Inf),
##'               M0 = list(ini = 0,
##'                       fixed = 1,
##'                       lower = 0,
##'                       upper = Inf),
##'                     residue =c( 0.0,  0.0,  0.0,  1.6,  4.0, 12.3, 13.5,
##' 12.7, 11.4, 11.6, 10.9,  9.5,  7.6))                           )
##' fit1 <- mkinfit.gui(SFO_SFO_gui)
##' summary(fit1)
##' fit2 <- IRLSkinfit.gui(SFO_SFO_gui)
##' summary(fit2)
##' }
##'
summary.kingui0 <- function(object, data = TRUE, distimes = TRUE, ff = TRUE, cov = FALSE,version="1.2011.701.11155",...) {
    options(warn=-1)
    param  <- object$par
  pnames <- names(param)
  p      <- length(param)
  covar  <- object$covar
  #rownames(covar) <- colnames(covar) <-pnames
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  se     <- sqrt(diag(covar) * resvar)
  #browser()
  lci <- param-qnorm(0.975)*se
  uci <- param+qnorm(0.975)*se
  names(se) <- pnames
  tval      <- param / se
  modVariance <- object$ssr / length(object$residuals)

  if (!all(object$start$lower >=0)) {
    message <- "Note that the one-sided t-test may not be appropriate if
      parameter values below zero are possible."
    warning(message)
  } else message <- "ok"

  param <- cbind(param, se, lci,uci, pt(tval, rdf, lower.tail = FALSE))
 ## adding one line when there is fixed by KinGui object
    pnames1 <- pnames

    fid <- which(object$fixed$by=='KinGui')
    if(length(fid)>0){
        nfid <- length(fid)
        param <- rbind(param,(cbind(object$fixed$value[fid],matrix(NA,nfid,4))))
        pnames1 <- c(pnames1, rownames(object$fixed)[fid])
    }
    dimnames(param) <- list(pnames1, c("Estimate", "Std. Error",'Lower CI','Upper CI',
                                      "Pr(>t)"))
  ##dimnames(param) <- list(pnames, c("Estimate", "Std. Error",'Lower CI','Upper CI',
    ##                                 "Pr(>t)"))
  if(cov)
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf), cov.unscaled = covar,
                cov.scaled = covar * resvar,
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param,version=version,pnames=pnames)
  else
    ans <- list(residuals = object$residuals,
                residualVariance = resvar,
                sigma = sqrt(resvar),
                modVariance = modVariance,
                df = c(p, rdf),
                info = object$info, niter = object$iterations,
                stopmess = message,
                par = param,version=version,pnames=pnames)

  ans$diffs <- object$diffs
  if(data) ans$data <- object$data
ans$start0 <- object$start0
ans$fixed0 <- object$fixed0
  ans$start <- object$start
  ans$fixed <- object$fixed
  ans$errmin <- object$errmin
  if(distimes) ans$distimes <- object$distimes
  if(ff) ans$ff <- object$ff
  class(ans) <- c('summary.kingui0', "summary.modFit")
  options(warn=0)
  return(ans)
}

# Expanded from print.summary.modFit
##' Print method for \code{\link{summary.kingui}}
##'
##' ##' Expanded from \code{\link{print.summary.modFit}} and
##' \code{\link{print.summary.mkinfit}}
##' @param x An object of class \code{summary.kingui}
##' @param digits How many digits should be printed after the decimal point.
##' @param detailed Ignore for now.
##' @param ... Other parameters to be passed into \code{\link{summary.kingui}}
##' @return \code{NULL}
##' @author Zhenglei Gao
##' @keywords internal
print.summary.kingui0 <- function(x, digits = max(3, getOption("digits") - 3),detailed=FALSE, ...) {
    cat(paste('Version:',x$version,'\n'))
    ## cat("\nVersion: ")
    ## print(x$version)
    ##  cat("\n")
    xx <- x[["diffs"]]
  cat("\nEquations:\n")
  for(i in 1:length(xx)) print(noquote(as.character(xx[i])))
  df  <- x$df
  rdf <- df[2]

  cat("\nStarting values for optimised parameters:\n")
    print(x$start0)

  cat("\nFixed parameter values:\n")

        if(length(x$fixed0$value) == 0) cat("None\n")
        else print(x$fixed0)


  cat("\nOptimised parameters:\n")
  printCoefmat(x$par, digits = digits, ...)

  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  printff <- !is.null(x$ff)
  if(printff){
    cat("\nEstimated formation fractions:\n")
    print(data.frame(ff = x$ff), digits=digits,...)
  }




  cat("\nChi2 error levels in percent :\n")
  print(x$errmin[,1:3], digits=digits,...)

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(myformat(x$distimes, digits=digits,scientific=F,...))
  }





  printcor <- (!is.null(x$cov.unscaled))
  if (printcor){
      cat("\nAdditional Statistics:\n")
      print(x$errmin[,4:ncol(x$errmin)], digits=digits,...)
      Corr <- try(cov2cor(x$cov.unscaled),silent=TRUE)
      if(!is.numeric(Corr))
      {
          warning('Covariance matrix cannot be calculated.')
          np <- nrow(x$cov.unscaled)
         Corr <- matrix(data = NA, nrow = np, ncol = np)
      }
      rownames(Corr) <- colnames(Corr) <- x$pnames ##rownames(x$par)
      cat("\nParameter correlation:\n")
      print(Corr, digits = digits, ...)
  }


  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(myformat(x$data, digits = digits, scientific=F,...), row.names = FALSE)
    #sprintf("%.4f",x$data)
  }

  invisible(x)
}
########################################
##'This function calculates the smallest relative error resulting
##' in passing the chi-squared test as defined in the FOCUS kinetics
##' report from 2006.
##'
##' This function is used internally by \code{\link{mkinfit.gui}},
##' \code{\link{IRLSkinfit.gui}}, and  \code{\link{mcmckinfit.gui}}.
##' @param errdata A data frame
##' @param n.parms Number of parameters
##' @param alpha Significance level
##' @keywords internal
##' @note Deprecated!
chi2err.bug <- function(errdata, n.parms, alpha = 0.05)
{
    means.mean <- mean(errdata$value_mean, na.rm = TRUE)
    df = length(errdata$value_mean) -sum(is.na(errdata$value_mean)) - n.parms
    C <- errdata$value_pred
    O <- errdata$value_mean
    b <-  sum((C-O)^2,na.rm=TRUE)
    a <- b/(means.mean)^2
    #a <- sum((C-O)^2/(means.mean)^2,na.rm=TRUE)
    err.sig <-sqrt(b/qchisq(1-alpha,df))
    err.min <- sqrt(1/qchisq(1-alpha,df)*a)
    O1 <- errdata$value_obs
    ef <- EF(O1,C)
    RMSE <- rmse(O1,C)
    R2 <- cor(O1,C,use='na.or.complete')^2
    return(list(err.min = 100*err.min, n.optim = n.parms, df = df,err.sig=err.sig,RMSE=RMSE,EF=ef,R2=R2))
}
