cost1<-function (P, state.ini.optim, state.ini.fixed, state.ini.optim.boxnames, 
    state.ini.fixed.boxnames, parms.fixed, outtimes, mkindiff, 
    mkinmod, observed, err, weight, scaleVar, pnames) 
{
    names(P) <- pnames
    if (length(state.ini.optim) > 0) {
        odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
        names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    }
    else {
        odeini <- state.ini.fixed
        names(odeini) <- c(state.ini.fixed.boxnames)
    }
    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], 
        parms.fixed)
    out <- ode(y = odeini, times = outtimes, func = mkindiff, 
        parms = odeparms)
    out_transformed <- data.frame(time = out[, "time"])
    for (var in names(mkinmod$map)) {
        if (length(mkinmod$map[[var]]) == 1) {
            out_transformed[var] <- out[, var]
        }
        else {
            out_transformed[var] <- rowSums(out[, mkinmod$map[[var]]])
        }
    }
    assign("out_predicted", out_transformed, inherits = TRUE)
    mC <- modCost(out_transformed, observed, y = "value", err = err, 
        weight = weight, scaleVar = scaleVar)
    return(mC)
}