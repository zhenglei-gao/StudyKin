plotobj <- function (mkinmodini,fit,eigen=FALSE,ctr=kingui.control(),plottitle='',
   err = NULL, weight = "none",
    scaleVar = FALSE,plotrange=c(0.9,1.1),plotn=100,plotparm=c(1,2),update=NULL,...)
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
    ##############

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

    return(mC)
}
    x0 <-fit$par[plotparm[1]]
    y0 <- fit$par[plotparm[2]]
    x <- seq(plotrange[1]*x0,plotrange[2]*x0,length=plotn)
    x <- rep(x,each=plotn)
    y <- seq(plotrange[1]*y0,plotrange[2]*y0,length=plotn)
    y <- rep(y,plotn)
    #DF[rep(seq(nrow(DF)), each = 3), ]
    pmat <- matrix(rep(fit$par,plotn*plotn),plotn*plotn,length(fit$par),byrow=TRUE)
    pmat[,plotparm] <- cbind(x,y)
    #browser()
   # observed <-data.frame(time=fit$data$time,name=fit$data$variable,value=fit$data$observed,err=fit$data$'err std')
    if(!is.null( fit$data$'err std'))observed$err <- fit$data$'err std'
    if(!is.null( fit$data$'err-std'))observed$err <- fit$data$'err-std'

    z <- apply(pmat,1,function(x){
        names(x) <- c( names(state.ini.optim) , names(parms.optim) )
        f <- cost(x)
        f <- f$model
        })
    ###

    plot3d(x,y,z)
    return(data.frame(x,y,z))
}
