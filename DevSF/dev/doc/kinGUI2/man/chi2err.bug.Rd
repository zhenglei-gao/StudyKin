\name{chi2err.bug}
\alias{chi2err.bug}
\title{This function calculates the smallest relative error resulting
in passing the chi-squared test as defined in the FOCUS kinetics
report from 2006.}
\usage{
  chi2err.bug(errdata, n.parms, alpha = 0.05)
}
\arguments{
  \item{errdata}{A data frame}

  \item{n.parms}{Number of parameters}

  \item{alpha}{Significance level}
}
\description{
  This function is used internally by
  \code{\link{mkinfit.gui}}, \code{\link{IRLSkinfit.gui}},
  and \code{\link{mcmckinfit.gui}}.
}
\note{
  Deprecated!
}
\keyword{internal}

