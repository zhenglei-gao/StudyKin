\name{KINGUII-package}
\alias{KINGUII-package}
\alias{KINGUII}
\docType{package}
\title{
Routines for fitting kinetic models NLS, IRLS and MCMC methods.
}
\description{
R-package KINGUII is a GUI version of package \code{\link{mkin}}. It
contains functions to fit kinetic models using NLS,IRLS and MCMC
methods. It also includes a function to set up the kineticmodel with
differential equations.

The main functions are:
  \itemize{
    \item Function \code{mkinmod.gui} to set up the kinetic model, the data, and the
  initial values.
    \item Functions to allow fitting of the model to
    data: \code{mkinfit.gui}, \code{IRLSkinfit.gui}, and \code{mcmckinfit.gui}.
    \item Functions to summarize the results and produce output which
    can be read by the GUI: \code{summary.kingui},
   \code{summary.mcmckingui}, \code{kingraph}, etc.
    }
}
\details{
\tabular{ll}{
Package: \tab KINGUII\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-06-21\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Johannes Ranke and Zhenglei Gao

Maintainer: Johannes Ranke <jranke@users.berlios.de>
}
\references{

}

\keyword{ package }
\seealso{
 \code{\link[mkin]{<package-mkin>}}
}
\examples{
\dontrun{
SFO_SFO_gui <- mkinmod.gui(Parent = list(type = "SFO", to = "Metab", sink = TRUE,
                           k = list(ini = 0.1,
                      fixed = 0,
                      lower = 0,
                      upper = Inf),
              M0 = list(ini = 195,
                      fixed = 0,
                      lower = 0,
                      upper = Inf),
                           FF = list(ini = c(.1),
                      fixed = c(0),
                      lower = c(0),
                      upper = c(1)),
                      time=c(0.0,2.8,   6.2,  12.0,  29.2,  66.8,  99.8,
127.5, 154.4, 229.9, 272.3, 288.1, 322.9),
                    residue = c( 157.3, 206.3, 181.4, 223.0, 163.2,
144.7,  85.0,  76.5,  76.4,  51.5,  45.5,  47.3, 42.7)),
                           Metab = list(type = "SFO",
                           k = list(ini = 0.1   ,
                      fixed = 0,
                      lower = 0,
                      upper = Inf),
              M0 = list(ini = 0,
                      fixed = 1,
                      lower = 0,
                      upper = Inf),
                    residue =c( 0.0,  0.0,  0.0,  1.6,  4.0, 12.3, 13.5,
12.7, 11.4, 11.6, 10.9,  9.5,  7.6))                           )
fit1 <- mkinfit.gui(SFO_SFO_gui)
fit2 <- IRLSkinfit.gui(SFO_SFO_gui)
fit3 <- mcmckinfit.gui(SFO_SFO_gui)
}
}
