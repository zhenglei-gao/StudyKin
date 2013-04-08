##' Fit a kinetic model using the IRLS algorithm.
##'
##' Instead of implicitly assuming equal error variances or giving arbitrary weights decided by the researcher as in the NLS algorithm,  an iteratively reweighted least squares (IRLS) algorithm was implemented to obtain the maximum likelihood estimates of the kinetic model parameters.
##' @title Fit a kinetic model using the IRLS algorithm.
##' @param mkinmodini  A list of class \code{\link{mkinmod.full}}, containing the kinetic model to be fitted to the data, and the initial parameter values, the observed data.
##' @param eigen  If TRUE,  the solution of the system of differential equation should be based on the spectral decomposition of the coefficient matrix in cases that this is possible.
##' @param plot If TRUE,the observed values and the numerical solutions should be plotted at each stage of the optimisation.
##' @param plottitle The title of the plot for visualizing the optimization process.
##' @param quiet  If TRUE, suppress printing out the current model cost after each(>1) improvement.
##' @param err  See argumetns of \code{\link{mkinfit.full}}
##' @param weight See argumetns of \code{\link{mkinfit.full}}
##' @param scaleVar See argumetns of \code{\link{mkinfit.full}}
##' @param ctr A list of control values for the estimation algorithm to replace the default values including maximum iterations and absolute error tolerance.  Defaults to the output of \code{\link{kingui.control}}.
##' @param irls.control A list of control values for the estimation algorithm to replace the default values including the maximum number of iterations for the outer iteration and the error tolerance level for the error variance estimation updating.
##' @param update If not NULL, should be a list of starting values obtained from other optimization methods.
##' @param useHsolnp Whether to use the hessian matrix derived from the solnp optimization algorithm.
##' @param ...  Further arguments that will be passed to \code{\link{modFit}}.
##' @return A list with  "kingui", "mkinfit" and "modFit" in the class attribute. A summary can be obtained by \code{\link{summary.kingui}}.
##' @author Zhenglei Gao
##' @examples
##' complex <- mkinmod.full(
##'  parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
##'  A1 = list(type = "SFO", to = "A2"),
##'  B1 = list(type = "SFO"),
##'  C1 = list(type = "SFO"),
##'  A2 = list(type = "SFO"),
##'  inpartri='default',
##'  outpartri='default',
##'  data=schaefer07_complex_case,
##'  weight=NULL)
##' Fit    <- IRLSkinfit.full(
##'            complex,
##'               plot      = TRUE,
##'               quiet     = TRUE,
##'               ctr       = kingui.control
##'                            (method = 'solnp',
##'                            submethod = 'Port',
##'                              maxIter = 100,
##'                            tolerance = 1E-06,
##'                            odesolver = 'lsoda'),
##'            irls.control = list(
##'                              maxIter = 10,
##'                            tolerance = 0.001))
##' @keywords Kinetic-Evaluations
##' @export
KinEval <- function(mkinmodini,
                    evalMethod=c('NLS','IRLS','MCMC','directMLE'),
                    optimMethod=c("LM","port","nls2","Nash","Marq", "Port", "Newton", "Nelder-Mead", "BFGS", 
					"CG","L-BFGS-B", "SANN", "Pseudo",'trust',
					'spg','ucminf','nmk','Rcgmin','Rvmmin','deoptim','solnp'),
                    eigen = FALSE,
                    plot = FALSE, plottitle='',quiet = FALSE,
                    ctr=kingui.control(),irls.control=list(),
                    update=NULL,useHsolnp=FALSE,...)
{
  ## Must be simple and clear this time. :)
  
  ##################### 
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
  
  ## -----------------------------------
  ## Get the parametrization.
  inpartri <- mkinmodini$inpartri
  outpartri <- mkinmodini$outpartri
  ##
  
  ## mkinmodini is an object by mkinmod.full
  parms.ini <- mkinmodini$parms.ini
  state.ini <- mkinmodini$state.ini
  lower <- mkinmodini$lower
  upper <- mkinmodini$upper
  fixed_parms <- mkinmodini$fixed_parms
  fixed_initials <- mkinmodini$fixed_initials
  mod_vars <- names(mkinmodini$diffs)
  observed <-  mkin_wide_to_long(mkinmodini$residue,time='time')
  observed$err <-c(as.matrix(mkinmodini$weightmat))
  ## Subset dataframe with mapped (modelled) variables
  observed <- subset(observed, name %in% names(mkinmodini$map))
  ## Get names of observed variables
  ## NOTE HERE: the order may not be the same as the input mkinmod.full differential equations list. ## XXXXX TODO XXXX Reorder them maybe a good idea if the data is given from a data file while the mkinmod.full is defined not following the colnames order, although it is already taken care of in the cost(P) function to reorder the odeini using mod_vars
  obs_vars = unique(as.character(observed$name))
  
  
  ## Name the parameters if they are not named yet ## usually they are already names
  if(is.null(names(parms.ini))) names(parms.ini) <- mkinmodini$parms
  
  ## Name the inital parameter values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars
  
  ## Parameters to be optimised
  parms.fixed <- parms.ini[fixed_parms]
  optim_parms <- setdiff(names(parms.ini), fixed_parms)
  parms.optim <- parms.ini[optim_parms]
  
  
  ## # ### ### ### ### ###
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
  
  eigen <- FALSE
  oldparms <- c(state.ini.optim,parms.optim)
  if (length(mkinmodini$map) == 1) {
    solution = "analytical"
  } else {
    if (is.matrix(mkinmodini$coefmat) & eigen) solution = "eigen"
    else solution = "deSolve"
  }
  ## always define mkindiff function since most of the time we will use it.
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
  
  
  # -----------------------------------
  ## Get the evaluation and optimization method.
  evalMethod <- match.arg(evalMethod)
  optimMethod <- match.arg(optimMethod)
  runIRLS <- TRUE ## determin whther to run an IRLS step especially before the MCMC step
  ## if the parameter values are provided to the MCMC, then do not run 
  if(evalMethod == 'MCMC' & (!is.null(update))) runIRLS <- FALSE  
  if(evalMethod == "NLS") runIRLS <- FALSE
  # ----------------------------------------
  # options(warn=-1) ## turn off the warning
  # ----------------------------------------  
  # -----------------------------------------
  ## Run the Optimization ##
  ## Always run a first step optimization.
  environment(kin_mod) <- environment()
  oldparms <- c(state.ini.optim,parms.optim)
  y<-observed$value
  if(optimMethod=="port") res0 <- nls(y ~ kin_mod(P,pnames=names(oldparms)),start=list(P=oldparms),lower=lower,upper=upper,algorithm="port")
  f <- function(P){
    observed$value-kin_mod(P)
  }
  if(optimMethod=="LM") res0 <- nls.lm(par=oldparms,lower=lower,upper=upper,fn=f)
  if(optimMethod== "Nash") {
    ## res0 <- nlxb(y ~ kin_mod(P,pnames=names(oldparms)),start=list(P=oldparms),lower=lower,upper=upper)
    res0 <- nlfb(start=oldparms, resfn=f,lower=lower,upper=upper)
  }
  if(optimMethod == "nls2"){
    P <- oldparms
    res0 <- nls2(y ~ kin_mod(P,pnames=names(oldparms)),start=t(oldparms),lower=lower,upper=upper)
  }
  ## if(evalMethod='MCMC'), stop here.
  if(evalMethod=='IRLS' & runIRLS ){ ## when IRLS, need an iterative step:
  
  }
  if(evalMethod=='MCMC'){ ## when MCMC, always use the different 
  ## 'sigma' setting and start with the fitted parameters from IRLS results
  }
  # ----------------------------------------
  
  # -----------------------------------------
  ## Result Summary Section ##
  # -----------------------------------------
  return(res0)
}


