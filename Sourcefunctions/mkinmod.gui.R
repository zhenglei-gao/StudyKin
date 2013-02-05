# $Id: mkinmod.gui.R

# Copyright (C) 2010 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>



### Example Usage




mkinmod.gui <- function(...)
{
  spec <- list(...)
  obs_vars <- names(spec)

  # The returned model will be a list of vectors, containing
  # differential equations, parameter names and a mapping from model variables
  # to observed variables, and all the parameter/state initial values, fixed or not, lower and upper bounds . If possible, a matrix representation of the
  # differential equations is included.

  parms <- vector()
  parms.ini <- vector()
  parms.lower <- vector()
  parms.upper <- vector()
  state.ini <- vector()
  state.lower <- vector()
  state.upper <- vector()
  diffs <- vector()
  ff <- vector()
  residue <- NULL
  weightmat <- NULL
  map <- list()
  fixed_parms <- NULL
  fixed_flag <- NULL
  fixed_initials <- NULL
  lower <- NULL
  upper <- NULL
  if(spec[[1]]$type %in% c("FOMC", "DFOP", "HS")) {
    mat = FALSE
  } else mat = TRUE
  time <- spec[[1]]$time
  residue <- time
  ############### Define the coefmat matrix/Jacobian of the ODE system
   if (mat) {
       n <- length(obs_vars)
       m <- matrix(nrow=n, ncol=n, dimnames=list(obs_vars, obs_vars))
  }



  # Establish list of differential equations
  for (varname in obs_vars)
  {
    if(is.null(spec[[varname]]$type)) stop(
      "Every argument to mkinmod must be a list containing a type component")
    if(!spec[[varname]]$type %in% c("SFO", "FOMC", "DFOP", "HS", "SFORB")) stop(
      "Available types are SFO, FOMC, DFOP, HS and SFORB only")
    new_parms <- vector()
    new_parms.ini <- vector()
    new_parms.lower <- vector()
    new_parms.upper <- vector()
    new_ff <- vector()
    if(is.null(spec[[varname]]$M0)) state.ini <- c(state.ini,0) else state.ini <- c(state.ini,spec[[varname]]$M0$ini)
    if(is.null(spec[[varname]]$M0)) state.lower <- c(state.lower,0) else state.lower <- c(state.lower,spec[[varname]]$M0$lower)
    if(is.null(spec[[varname]]$M0)) state.upper <- c(state.upper,Inf) else state.upper <- c( state.upper,spec[[varname]]$M0$upper)
    if(spec[[varname]]$M0$fixed==1 || is.null(spec[[varname]]$M0)) fixed_initials <- c(fixed_initials,varname)
    residue <- cbind(residue,spec[[varname]]$residue)
    if(!is.null(spec[[varname]]$weight)) weightmat <- cbind(weightmat,spec[[varname]]$weight) else weightmat <- cbind(weightmat,rep(1,length(spec[[varname]]$residue)))
    # New (sub)compartments (boxes) needed for the model type
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      FOMC = varname,
      DFOP = varname,
      HS = varname,
      SFORB = paste(varname, c("free", "bound"), sep="_")
    )
    map[[varname]] <- new_boxes
    names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes))

    # Start a new differential equation for each new box
    new_diffs <- paste("d_", new_boxes, " =", sep="")

    # Turn on sink if not specified otherwise
    if(is.null(spec[[varname]]$sink)) spec[[varname]]$sink <- TRUE

    # Construct and add FOMC term and add FOMC parameters if needed
    if(spec[[varname]]$type == "FOMC") {
      if(match(varname, obs_vars) != 1) {
        stop("Type FOMC is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
         warning("Turning off the sink for the FOMC model is not fully tested")
        }
      # From p. 53 of the FOCUS kinetics report
      alphaname<-paste("alpha", new_boxes[[1]],  sep="_")
      betaname<-paste("beta", new_boxes[[1]],  sep="_")

      ##nonlinear_term <- paste("(alpha/beta) * ((time/beta) + 1)^-1 *", new_boxes[[1]])
      nonlinear_term <- paste("(",alphaname,"/",betaname,") * ((time/",betaname,") + 1)^-1 *", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      ##new_parms <- c("alpha", "beta")
      new_parms <- c(alphaname, betaname)
      new_parms.ini <- c(spec[[varname]]$alpha$ini,spec[[varname]]$beta$ini)
      new_parms.lower <- c(spec[[varname]]$alpha$lower,spec[[varname]]$beta$lower)
      new_parms.upper <- c(spec[[varname]]$alpha$upper,spec[[varname]]$beta$upper)
      ##if(spec[[varname]]$alpha$fixed==1) fixed_parms <- c(fixed_parms,'alpha')
      ##if(spec[[varname]]$beta$fixed==1) fixed_parms <- c(fixed_parms,'beta')
      if(spec[[varname]]$alpha$fixed==1) fixed_parms <- c(fixed_parms,alphaname)
      if(spec[[varname]]$beta$fixed==1) fixed_parms <- c(fixed_parms,betaname)
      #new_ff <- spec[[varname]]$FF$ini
    }

    # Construct and add DFOP term and add DFOP parameters if needed
    if(spec[[varname]]$type == "DFOP") {
      if(match(varname, obs_vars) != 1) {
        stop("Type DFOP is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
        warning("Turning off the sink for the DFOP model is not fully tested!")
      }
      # From p. 57 of the FOCUS kinetics report
      #nonlinear_term <- paste("((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) *", new_boxes[[1]])
      k1name <-paste("k1", new_boxes[[1]],  sep="_")
      k2name <- paste("k2", new_boxes[[1]],  sep="_")
      gname <- paste("g", new_boxes[[1]],  sep="_")
      nonlinear_term <- paste("((",k1name,"*", gname ,"*", "exp(-",k1name, "* time) + ",k2name, "* (1 - ",gname,") * exp(-",k2name, "* time)) / (",gname, "* exp(-",k1name,"* time) + (1 - ",gname,") * exp(-",k2name, "* time))) *", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      #new_parms <- c("k1", "k2", "g")
      new_parms <- c(k1name,k2name,gname)
      new_parms.ini <- c(spec[[varname]]$k1$ini,spec[[varname]]$k2$ini,spec[[varname]]$g$ini)
      new_parms.lower <- c(spec[[varname]]$k1$lower,spec[[varname]]$k2$lower,spec[[varname]]$g$lower)
      new_parms.upper <- c(spec[[varname]]$k1$upper,spec[[varname]]$k2$upper,spec[[varname]]$g$upper)
      ## if(spec[[varname]]$k1$fixed==1) fixed_parms <- c(fixed_parms,'k1')
      ## if(spec[[varname]]$k2$fixed==1) fixed_parms <- c(fixed_parms,'k2')
      ## if(spec[[varname]]$g$fixed==1) fixed_parms <- c(fixed_parms,'g')
      if(spec[[varname]]$k1$fixed==1) fixed_parms <- c(fixed_parms,k1name)
      if(spec[[varname]]$k2$fixed==1) fixed_parms <- c(fixed_parms,k2name)
      if(spec[[varname]]$g$fixed==1) fixed_parms <- c(fixed_parms,gname)

      #new_ff <- spec[[varname]]$FF$ini
    }

    # Construct and add HS term and add HS parameters if needed
    if(spec[[varname]]$type == "HS") {
      if(match(varname, obs_vars) != 1) {
        stop("Type HS is only possible for the first compartment, which is assumed to be the source compartment")
      }
      if(spec[[varname]]$sink == FALSE) {
        warning("Turning off the sink for the HS model is not fully tested")
      }
      # From p. 55 of the FOCUS kinetics report
      k1name<-paste("k1", new_boxes[[1]],  sep="_")
      k2name<-paste("k2", new_boxes[[1]],  sep="_")
      tbname<-paste("tb", new_boxes[[1]],  sep="_")

      #nonlinear_term <- paste("ifelse(time <= tb, k1, k2)", "*", new_boxes[[1]])
      nonlinear_term <- paste("ifelse(time <=",tbname,",",k1name,",", k2name,")", "*", new_boxes[[1]])
      new_diffs[[1]] <- paste(new_diffs[[1]], "-", nonlinear_term)
      #new_parms <- c("k1", "k2", "tb")
      new_parms <- c(k1name, k2name, tbname)
      new_parms.ini <- c(spec[[varname]]$k1$ini,spec[[varname]]$k2$ini,spec[[varname]]$tb$ini)
      new_parms.lower <- c(spec[[varname]]$k1$lower,spec[[varname]]$k2$lower,spec[[varname]]$tb$lower)
      new_parms.upper <- c(spec[[varname]]$k1$upper,spec[[varname]]$k2$upper,spec[[varname]]$tb$upper)
      ## if(spec[[varname]]$k1$fixed==1) fixed_parms <- c(fixed_parms,'k1')
      ## if(spec[[varname]]$k2$fixed==1) fixed_parms <- c(fixed_parms,'k2')
      ## if(spec[[varname]]$tb$fixed==1) fixed_parms <- c(fixed_parms,'tb')
      if(spec[[varname]]$k1$fixed==1) fixed_parms <- c(fixed_parms,k1name)
      if(spec[[varname]]$k2$fixed==1) fixed_parms <- c(fixed_parms,k2name)
      if(spec[[varname]]$tb$fixed==1) fixed_parms <- c(fixed_parms,tbname)
      #new_ff <- spec[[varname]]$FF$ini
    }
    if(spec[[varname]]$type %in% c("SFO")) {
        k_compound_sink <- paste("k", new_boxes[[1]],  sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        new_diffs[[1]] <- paste(new_diffs[[1]], sink_term)
        new_parms <- k_compound_sink
        #if(is.null(spec[[varname]]$FF)) new_parms.ini <-spec[[varname]]$k$ini else new_parms.ini <-spec[[varname]]$k$ini*(1-sum(spec[[varname]]$FF$ini))
        new_parms.ini <-spec[[varname]]$k$ini
        new_parms.lower <- spec[[varname]]$k$lower
        new_parms.upper <- spec[[varname]]$k$upper
        #if(spec[[varname]]$k$fixed==1 & is.null(spec[[varname]]$FF)) fixed_parms <- c(fixed_parms,k_compound_sink)## there is no need, prod(NULL)=1
        #if(spec[[varname]]$k$fixed==1 & prod(spec[[varname]]$FF$fixed)==1) fixed_parms <- c(fixed_parms,k_compound_sink)
        if(spec[[varname]]$k$fixed==1) fixed_parms <- c(fixed_parms,k_compound_sink)

        ##

    }

    # Construct terms for transfer to sink and add if appropriate

    if(spec[[varname]]$sink) {
      # Add first-order sink term to first (or only) box for SFO and SFORB

     # Add first-order sink term to first (or only) box for SFO and SFORB
      if(spec[[varname]]$type %in% c("SFORB")) {
        k_compound_sink <- paste("k", new_boxes[[1]], sep="_")
        sink_term <- paste("-", k_compound_sink, "*", new_boxes[[1]])
        new_diffs[[1]] <- paste(new_diffs[[1]], sink_term)
        new_parms <- k_compound_sink
        if(is.null(spec[[varname]]$FF)) new_parms.ini <-spec[[varname]]$k$ini else new_parms.ini <-spec[[varname]]$k$ini*(1-sum(spec[[varname]]$FF$ini))
        new_parms.lower <- spec[[varname]]$k$lower
        new_parms.upper <- spec[[varname]]$k$upper
        if(spec[[varname]]$k$fixed==1 & is.null(spec[[varname]]$FF)) fixed_parms <- c(fixed_parms,k_compound_sink)## there is no need, prod(NULL)=1
        if(spec[[varname]]$k$fixed==1 & prod(spec[[varname]]$FF$fixed)==1) fixed_parms <- c(fixed_parms,k_compound_sink)

        ##
      }


    # Add reversible binding if appropriate
    if(spec[[varname]]$type == "SFORB") {
      k_free_bound <- paste("k", varname, "free", "bound", sep="_")
      k_bound_free <- paste("k", varname, "bound", "free", sep="_")
      reversible_binding_terms <- c(
        paste("-", k_free_bound, "*", new_boxes[[1]], "+", k_bound_free, "*", new_boxes[[2]]),
        paste("+", k_free_bound, "*", new_boxes[[1]], "-", k_bound_free, "*", new_boxes[[2]]))
      new_diffs <- paste(new_diffs, reversible_binding_terms)
      new_parms <- c(new_parms, k_free_bound, k_bound_free)
      new_parms.ini <-c(new_parms.ini,spec[[varname]]$k_free_bound$ini,spec[[varname]]$k_bound_free$ini)
      new_parms.lower <-c(new_parms.lower,spec[[varname]]$k_free_bound$lower,spec[[varname]]$k_bound_free$lower)
      new_parms.upper <-c(new_parms.upper,spec[[varname]]$k_free_bound$lower,spec[[varname]]$k_bound_free$upper)
      if(spec[[varname]]$k_free_bound$fixed==1) fixed_parms <- c(fixed_parms,k_free_bound)
      if(spec[[varname]]$k_bound_free$fixed==1) fixed_parms <- c(fixed_parms,k_bound_free)


    }
      if(mat) m[varname,varname] <- paste("-k", varname, sep="_")
  }else{
      if(spec[[varname]]$type %in% c("SFO") && mat) m[varname,varname] <- paste("-k", varname,  sep="_")
  }

    # Add observed variable to model
    parms <- c(parms, new_parms)
    parms.ini <- c(parms.ini, new_parms.ini)
    parms.lower <- c(parms.lower, new_parms.lower)
    parms.upper <- c(parms.upper, new_parms.upper)
    #FF_value<- c(FF_value,new_ff)
    names(new_diffs) <- new_boxes
    diffs <- c(diffs, new_diffs)
}
   ## Add a start component for later report usage.
    start <- NULL

    ## if(outpartri=='default'){
        start_parms <- setdiff(parms, fixed_parms)
        names(parms.ini) <- parms
        names(parms.lower) <- parms
        names(parms.upper) <- parms
        start_fixed <- rep(0,length(parms))
        names(start_fixed) <- parms
        start_fixed[fixed_parms] <- 1
        start <- data.frame(initial=parms.ini,lower=parms.lower,upper=parms.upper,fixed=start_fixed)
        ##start <-data.frame(initial=parms.ini[start_parms],lower=parms.lower[start_parms],upper=parms.upper[start_parms])
        ##rownames(start) <- start_parms
    #################################################
  fixed_flag <- rep('user',length(fixed_parms))
  # Transfer between compartments
  for (varname in obs_vars) {
    to <- spec[[varname]]$to
    FF <- spec[[varname]]$FF$ini  ## formation fraction for every compartment
 if(spec[[varname]]$sink==FALSE)
         {
             ### Then the ff to the last compartment must be calculated and not as an optimization parameter.
             if(is.null(to))
             {
                 ## # In this case, the sink compartment cannot be turned off.
                 for(k in obs_vars)
                 {
                     if(k!=varname)
                     {
                         if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
                     }
                 }

             }
           if(!is.null(to)) {
               ## set 0 in the coefmat
               for(k in obs_vars)
               {
                   if(!(k %in% to))
                   {
                       if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
                   }
               }
               f <- ForwardCalcFF(FF) ## the starting values that will be used in the optimization program so that the optimization does not need to deal with the sum ff==1 contraint.###
        names(f) <-to
        origin_box <- switch(spec[[varname]]$type,
        SFO = varname,
        FOMC = varname,
        DFOP = varname,
        HS = varname,
        SFORB = paste(varname, "free", sep="_"))
      fraction_left <- NULL
               nto <- length(to)
      for (target in to) {
        index <- match(target,to) ### find formation fraction to a corresponding compartment
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
## ############################################################
         ## add starting values for the optimization with formation fractions
                    start <-rbind(start,c(FF[index],spec[[varname]]$FF$lower[index],spec[[varname]]$FF$upper[index],spec[[varname]]$FF$fixed[index]))
                    rownames(start)[nrow(start)] <- paste('ff',varname,'to',target,sep='_')
# ############################################################

        if(spec[[varname]]$type %in% c("SFORB")) {
          k_from_to <- paste("k", origin_box, target_box, sep="_")
          diffs[[origin_box]] <- paste(diffs[[origin_box]], "-",
            k_from_to, "*", origin_box)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+",
            k_from_to, "*", origin_box)
          parms <- c(parms, k_from_to)
          parms.ini<- c(parms.ini, spec[[varname]]$k$ini*FF[index])
          parms.lower <- c(parms.lower, spec[[varname]]$k$lower)
          parms.upper <- c(parms.upper, spec[[varname]]$k$upper)
          if(spec[[varname]]$k$fixed==1 & spec[[varname]]$FF$fixed[index]==1) {
              fixed_parms <- c(fixed_parms,k_from_to)
              fixed_flag <- c(fixed_flag,'user')
          }
        }
        if(spec[[varname]]$type %in% c("SFO")) {
            fraction_to_target = paste('f',origin_box,'to', target, sep="_")
            fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
            sep="")
            if(is.null(fraction_left)) {
                fraction_really_to_target = fraction_to_target
                fraction_left = fraction_not_to_target
            } else {
                fraction_really_to_target = paste(fraction_left, " * ",
                fraction_to_target, sep="")
                fraction_left = paste(fraction_left, " * ",
                fraction_not_to_target, sep="")
            }
            ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
            diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                  paste("k", origin_box, sep="_"),'*',ff[paste(origin_box,'to', target, sep="_")], "*", origin_box)
            parms <- c(parms, fraction_to_target)
            parms.ini <- c(parms.ini,f[target])
            if(spec[[varname]]$FF$fixed[index]==1) {
                fixed_parms <- c(fixed_parms,fraction_to_target)
                fixed_flag <- c(fixed_flag,'user')
            }
            ## IN NO SINK, then fixed parms should be adding 1, and the last tranformed formation fraction should be fixed at 1!!!!!!!!!!!!
            #browser()
            if(index==nto)
                {if(spec[[varname]]$FF$fixed[index]==0)
                 {
                    fixed_parms <- c(fixed_parms,fraction_to_target)
                     fixed_flag <- c(fixed_flag,'KinGui')
                }else{
                     warning('You need to switch the order if the formation fraction for the last to compartment is fixed at a certain value but you turn off the sink compartment and you have multiple to compartments!')
                     fixed_flag[length(fixed_flag)] <- 'KinGui'
                }
                    parms.ini[length(parms.ini)] <-1
                }

            ###############################################################
            parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
            parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
            if(mat){
                if(is.na(m[target,varname])) m[target,varname] <-paste(paste("k", origin_box, sep="_"),'*',fraction_really_to_target) else m[target,varname] <-paste(m[target,varname],'+',paste("k", origin_box, sep="_"),'*',fraction_really_to_target,sep='')
            }
        }

        if(spec[[varname]]$type %in% c("FOMC", "DFOP", "HS")) {

          fraction_to_target = paste('f',origin_box,'to', target, sep="_")
          fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
            sep="")
          if(is.null(fraction_left)) {
            fraction_really_to_target = fraction_to_target
            fraction_left = fraction_not_to_target
          } else {
            fraction_really_to_target = paste(fraction_left, " * ",
              fraction_to_target, sep="")
            fraction_left = paste(fraction_left, " * ",
              fraction_not_to_target, sep="")
          }
          ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
          #ff[target_box] = fraction_really_to_target
          #diffs[[target_box]] <- paste(diffs[[target_box]], "+",
          #  ff[target_box], "*", nonlinear_term)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+",
            fraction_really_to_target, "*", nonlinear_term)
          parms <- c(parms, fraction_to_target)
          parms.ini <- c(parms.ini, f[index])
          if(spec[[varname]]$FF$fixed[index]==1){
              fixed_parms <- c(fixed_parms,fraction_to_target)
              fixed_flag <- c(fixed_flag,'user')
          }
           ## IN NO SINK, then fixed parms should be adding 1!!!!!!!!!!!!
            #browser()
            if(index==nto)
                {if(spec[[varname]]$FF$fixed[index]==0)
                 {
                    fixed_parms <- c(fixed_parms,fraction_to_target)
                    fixed_flag <- c(fixed_flag,'KinGui')
                }else{
                    warning('You need to switch the order if the formation fraction for the last to compartment is fixed at a certain value but you turn off the sink compartment and you have multiple to compartments!')
                    fixed_flag[length(fixed_flag)] <- 'KinGui'

                }
                    parms.ini[length(parms.ini)] <-1
                }
          parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
          parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
        }
    }## end for(target in to)
           }## end if(!is.null(to))


         }else{
    if(is.null(to))
        {
            for(k in obs_vars)
            {
                if(k!=varname)
                {
                    if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
                }
            }

        }
    if(!is.null(to)) {
         for(k in obs_vars)
            {
                if(!(k %in% to))
                {
                    if(mat){if(is.na(m[k,varname])) m[k,varname] <- 0}
                }
            }
        f <- ForwardCalcFF(FF) ## the starting values that will be used in the optimization program o that the optimization does not need to deal with the sum ff==1 contraint.
        names(f) <-to
        origin_box <- switch(spec[[varname]]$type,
        SFO = varname,
        FOMC = varname,
        DFOP = varname,
        HS = varname,
        SFORB = paste(varname, "free", sep="_"))
      fraction_left <- NULL
      for (target in to) {
        index <- match(target,to) ### find formation fraction to a corresponding compartment
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          SFORB = paste(target, "free", sep="_"))
        if(spec[[varname]]$type %in% c("SFORB")) {
          k_from_to <- paste("k", origin_box, target_box, sep="_")
          diffs[[origin_box]] <- paste(diffs[[origin_box]], "-",
            k_from_to, "*", origin_box)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+",
            k_from_to, "*", origin_box)
          parms <- c(parms, k_from_to)
          parms.ini<- c(parms.ini, spec[[varname]]$k$ini*FF[index])
          parms.lower <- c(parms.lower, spec[[varname]]$k$lower)
          parms.upper <- c(parms.upper, spec[[varname]]$k$upper)
          if(spec[[varname]]$k$fixed==1 & spec[[varname]]$FF$fixed[index]==1) {
              fixed_parms <- c(fixed_parms,k_from_to)
              fixed_flag <- c(fixed_flag,'user')
          }
        }
        if(spec[[varname]]$type %in% c("SFO")) {
            fraction_to_target = paste('f',origin_box,'to', target, sep="_")
            fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
            sep="")
            if(is.null(fraction_left)) {
                fraction_really_to_target = fraction_to_target
                fraction_left = fraction_not_to_target
            } else {
                fraction_really_to_target = paste(fraction_left, " * ",
                fraction_to_target, sep="")
                fraction_left = paste(fraction_left, " * ",
                fraction_not_to_target, sep="")
            }
            ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
            diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                  paste("k", origin_box, sep="_"),'*',ff[paste(origin_box,'to', target, sep="_")], "*", origin_box)
            parms <- c(parms, fraction_to_target)
            parms.ini <- c(parms.ini,f[target])
            if(spec[[varname]]$FF$fixed[index]==1) {
                fixed_parms <- c(fixed_parms,fraction_to_target)
                fixed_flag <- c(fixed_flag,'user')
            }
            parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
            parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
            if(mat){
                if(is.na(m[target,varname])) m[target,varname] <-paste(paste("k", origin_box, sep="_"),'*',fraction_really_to_target) else m[target,varname] <-paste(m[target,varname],'+',paste("k", origin_box, sep="_"),'*',fraction_really_to_target,sep='')
            }
        }

        if(spec[[varname]]$type %in% c("FOMC", "DFOP", "HS")) {

          fraction_to_target = paste('f',origin_box,'to', target, sep="_")
          fraction_not_to_target = paste("(1 - ", fraction_to_target, ")",
            sep="")
          if(is.null(fraction_left)) {
            fraction_really_to_target = fraction_to_target
            fraction_left = fraction_not_to_target
          } else {
            fraction_really_to_target = paste(fraction_left, " * ",
              fraction_to_target, sep="")
            fraction_left = paste(fraction_left, " * ",
              fraction_not_to_target, sep="")
          }
          ff[paste(origin_box,'to', target, sep="_")] = fraction_really_to_target
          #ff[target_box] = fraction_really_to_target
          #diffs[[target_box]] <- paste(diffs[[target_box]], "+",
          #  ff[target_box], "*", nonlinear_term)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+",
            fraction_really_to_target, "*", nonlinear_term)
          parms <- c(parms, fraction_to_target)
          parms.ini <- c(parms.ini, f[index])
          if(spec[[varname]]$FF$fixed[index]==1) {
              fixed_parms <- c(fixed_parms,fraction_to_target)
              fixed_flag <- c(fixed_flag,'user')
          }
          parms.lower <- c(parms.lower, spec[[varname]]$FF$lower[index])
          parms.upper <- c(parms.upper, spec[[varname]]$FF$upper[index])
        }
		}## end for(target in to)
    }## end if(!is.null(to))
	}## end if(sink==FALSE){}else{}
}## end for(varname in ob_vars)
  names(parms.ini) <- parms
  names(parms.lower) <- parms
  names(parms.upper) <- parms
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]
  parms.lower <- parms.lower[optim_parms]
  parms.upper <- parms.upper[optim_parms]
  names(state.ini) <- obs_vars
  names(state.lower) <- obs_vars
  names(state.upper) <- obs_vars
  state.ini.fixed <- state.ini[fixed_initials]
  optim_initials <- setdiff(names(state.ini), fixed_initials)
  state.ini.optim <- state.ini[optim_initials]
  state.lower <- state.lower[optim_initials]
  state.upper <- state.upper[optim_initials]
  #fit <- modFit(cost, c(state.ini.optim, parms.optim), lower = lower, upper = upper, ...)
  lower <- c(state.lower,parms.lower)
  upper <- c(state.upper,parms.upper)
  colnames(residue) <- c('time',obs_vars)
  colnames(weightmat) <- obs_vars
   if(!is.null(start)) start$type <-  rep("deparm", nrow(start))
  model <- list(diffs = diffs, parms = parms, map = map,parms.ini=parms.ini,state.ini=state.ini,lower=lower,upper=upper,fixed_parms=fixed_parms,fixed_flag=fixed_flag,fixed_initials=fixed_initials,residue=as.data.frame(residue),weightmat=as.data.frame(weightmat),start=start)

  # Create coefficient matrix if appropriate
  if (mat) {

    model$coefmat <- m
  }

  if (exists("ff")) model$ff = ff
  class(model) <- "mkinmod"
  invisible(model)
}






###############################################################

BackCalcFF <- function(ff)
  {
    ### From the ff used in the program, back calculate the formation faction
    l <- length(ff)
    trueff <- ff
    if(l>1)
    {
        for(i in 1:l)
        {
            if(i>1)
            {
                for(j in 1:(i-1))
                    trueff[i] <- trueff[i]*(1-ff[j])
            }#else trueff[i] <- ff[i]
        }
    }
    trueff
  }

ForwardCalcFF <- function(trueff)
  {
    ### From the formation faction provided by the user(true formation fractions), calculte the f used in the program
    l <- length(trueff)
    ff <- trueff
    for(i in 1:l)
      {
        if(i>1)
          {
            for(j in 1:(i-1))
              ff[i] <- ff[i]/(1-ff[j])
          }#else ff[i] <- trueff[i]
      }
    ff
  }
