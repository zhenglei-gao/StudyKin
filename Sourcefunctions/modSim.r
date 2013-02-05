## a simulation function used for internal simluation only
## Example Usage:
## parms.ini <- c(9.163040e-02,    3.756059e-01,    8.669980e+00,
##4.090553e-18 ,1.766376e-02,    1.164488e-02,   3.314102e-01,    6.298495e-01,
##1.484640e-01,1.184215e-01,    1.729477e-05,    9.972716e-01,    2.134810e-02,    1.976447e-02 )
## state.ini<- c(100,rep(0,5))
## out1 <- modSim(MESO1,unique(observed2$time),1,parms.ini,state.ini)
nonegative <- function(x)
{
    x[x<0] <- 0
    x
}
modSim <- function(mkinmod,outtimes,sigma=1,parms.ini,state.ini,plot=TRUE,atol=1e-6,eigen=FALSE,plottitle='',fname=NULL,ghost=NULL) {

    mod_vars <- names(mkinmod$diffs)
    ## Name the parameters if they are not named yet
    if(is.null(names(parms.ini))) names(parms.ini) <- mkinmod$parms

    ## Name the inital parameter values if they are not named yet
    if(is.null(names(state.ini))) names(state.ini) <- mod_vars

    ## Decide if the solution of the model can be based on a simple analytical
    ## formula, the spectral decomposition of the matrix (fundamental system)
    ## or a numeric ode solver from the deSolve package
    if (length(mkinmod$map) == 1) {
        solution = "analytical"
    } else {
        if (is.matrix(mkinmod$coefmat) & eigen) solution = "eigen"
        else solution = "deSolve"
    }

    odeparms <- parms.ini
    odeini <- state.ini
    ## Create a function calculating the differentials specified by the model
    ## if necessary
    if(solution == "deSolve") {
        mkindiff <- function(t, state, parms) {
            time <- t
            diffs <- vector()
            for (box in mod_vars)
            {
                diffname <- paste("d", box, sep="_")
                diffs[diffname] <- with(as.list(c(time,state, parms)),
                                        eval(parse(text=mkinmod$diffs[[box]])))
            }
            return(list(c(diffs)))
        }
    }

    # Solve the system
    if (solution == "analytical") {
        parent.type = names(mkinmod$map[[1]])[1]
        parent.name = names(mkinmod$diffs)[[1]]
        o <- switch(parent.type,
                    SFO = SFO.solution(outtimes,
                    evalparse(parent.name),
                    evalparse(paste("k", parent.name, "sink", sep="_"))),
                    FOMC = FOMC.solution(outtimes,
                    evalparse(parent.name),
                    evalparse("alpha"), evalparse("beta")),
                    DFOP = DFOP.solution(outtimes,
                    evalparse(parent.name),
                    evalparse("k1"), evalparse("k2"),
                    evalparse("g")),
                    HS = HS.solution(outtimes,
                    evalparse(parent.name),
                    evalparse("k1"), evalparse("k2"),
                    evalparse("tb")),
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
        coefmat.num <- matrix(sapply(as.vector(mkinmod$coefmat), evalparse),
                              nrow = length(mod_vars))
        e <- eigen(coefmat.num)
        c <- solve(e$vectors, odeini)
        f.out <- function(t) {
            e$vectors %*% diag(exp(e$values * t), nrow=length(mod_vars)) %*% c
        }
        o <- matrix(mapply(f.out, outtimes),
                    nrow = length(mod_vars), ncol = length(outtimes))
        dimnames(o) <- list(mod_vars, outtimes)
        out <- cbind(time = outtimes, t(o))
    }
    if (solution == "deSolve")
    {
        out <- ode(
                   y = odeini,
                   times = outtimes,
                   func = mkindiff,
                   parms = odeparms,
                   atol = atol
                   )
    }
    out <- as.data.frame(out)
    out0 <- out
    n <- ncol(out)-1
    m <- nrow(out)
    if(length(sigma)<n) sigma <- c(sigma,rep(1e-4,n-length(sigma)))
    for(i in 1:n)
    {
        out[,i+1] <- nonegative(out[,i+1]+rnorm(m,0,sigma[i]))
    }

    if(!is.null(ghost))
    {
        for(names in ghost)
    {
        out[,names] <- rep(NA,m)
    }
    }
    out1 <- mkin_wide_to_long(out,time='time')
    if(plot==TRUE)
    {
        outtimes_plot = seq(min(outtimes), max(outtimes), length.out=100)
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
         ## o_plot <- switch(parent.type,
          ##   SFO = SFO.solution(outtimes_plot,
          ##       evalparse(parent.name),
          ##       evalparse(paste("k", parent.name, "sink", sep="_"))),
          ##   FOMC = FOMC.solution(outtimes_plot,
          ##       evalparse(parent.name),
          ##       evalparse("alpha"), evalparse("beta")),
          ##   DFOP = DFOP.solution(outtimes_plot,
          ##       evalparse(parent.name),
          ##       evalparse("k1"), evalparse("k2"),
          ##       evalparse("g")),
          ##   HS = HS.solution(outtimes_plot,
          ##       evalparse(parent.name),
          ##       evalparse("k1"), evalparse("k2"),
          ##       evalparse("tb")),
          ##   SFORB = SFORB.solution(outtimes_plot,
          ##       evalparse(parent.name),
          ##       evalparse(paste("k", parent.name, "bound", sep="_")),
          ##       evalparse(paste("k", sub("free", "bound", parent.name), "free", sep="_")),
          ##       evalparse(paste("k", parent.name, "sink", sep="_")))
          ## )
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
        for (var in names(mkinmod$map)) {
          if((length(mkinmod$map[[var]]) == 1) || solution == "analytical") {
            out_transformed_plot[var] <- out_plot[, var]
          } else {
            out_transformed_plot[var] <- rowSums(out_plot[, mkinmod$map[[var]]])
          }
        }
        obs_vars <- names(out)[-1]
        plot(0, type="n",
          xlim = range(out$time), ylim = range(out1$value, na.rm=TRUE),
          xlab = "Time", ylab = "Observed",main=plottitle)
        col_obs <- pch_obs <- 1:length(obs_vars)
        names(col_obs) <- names(pch_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(out1, name == obs_var, c(time, value)),
            pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_transformed_plot$time, out_transformed_plot[-1])
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars,
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
    }
    if(!is.null(fname)) write.table(out,fname,sep='\t')
    mod <- mkinmod
    mod$residue <-out
    mod$weightmat <- matrix(1,m,n)
    colnames(mod$weightmat) <- colnames(mkinmod$weightmat)
    mod$state.ini <- c(100,rep(0,n-1))
    names(mod$state.ini) <- names(mkinmod$state.ini)
    mod$parms.ini <-rep(0.1,length(parms.ini))
    names(mod$parms.ini) <- names(mkinmod$parms.ini)
    mod$parms.ini[mkinmod$fixed_parms] <-mkinmod$parms.ini[mkinmod$fixed_parms]
    mod$state.ini[mkinmod$fixed_initials] <-mkinmod$state.ini[mkinmod$fixed_initials]
    return(list(data=out,mod=mod,data0=out0))
}
