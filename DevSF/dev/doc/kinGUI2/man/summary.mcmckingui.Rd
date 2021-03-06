\name{summary.mcmckingui}
\alias{summary.mcmckingui}
\title{S3 summary method for class \code{mcmckingui}}
\usage{
  summary(object, remove = NULL,
    data=TRUE,distimes=TRUE,ff=TRUE,
    version="1.2011.922.1530",...)
}
\arguments{
  \item{object}{A fitted object of class \code{mcmckingui}
  from the result of \code{\link{mcmckinfit.full}}.}

  \item{remove}{The iterations should be removed from
  summary statistics calculations.}

  \item{data}{If TRUE, include in the returned values a
  data frame containing the observed and predicted values
  with residuals and estimated standard deviations or
  weights.}

  \item{distimes}{If TRUE, DT50 and DT90 values should be
  included.}

  \item{ff}{If TRUE, the formation fraction should be
  calculated from the estimated transformed parameters.}

  \item{version}{A version number indicating which version
  of the fit function has been used.}

  \item{...}{Optional arguments passed to methods like
  \code{print}.}
}
\value{
  The summary function returns a list of components from
  the results of the optimizat ion routine used.
}
\description{
  Lists model equations, the chi2 error levels calculated
  according to FOCUS guidance (2006) and optionally the
  data, consisting of observed, predicted and residual
  values, the correlation matrix and also the same summary
  statistics for class \code{modMCMC}
}
\author{
  Zhenglei Gao
}
\seealso{
  \code{\link{summary.kingui}}.
}

