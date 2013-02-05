\name{summary.kingui0}
\alias{summary.kingui0}
\title{S3 method for class 'kingui0'}
\usage{
  summary.kingui0(object, data = TRUE, distimes = TRUE,
    ff = TRUE, cov = FALSE, version = "1.2011.701.11155",
    ...)
}
\arguments{
  \item{object}{An object of class 'kingui' from the result
  of NLS or IRLS fit.}

  \item{data}{If TRUE, include in the returned values a
  data frame containing the observed and predicted values
  with residuals and estimated standard deviations or
  weights.}

  \item{distimes}{If TRUE, DT50 and DT90 values should be
  included.}

  \item{ff}{If TRUE, the formation fraction should be
  calculated from the estimated transformed parameters.}

  \item{cov}{If TRUE, parameter covariances should be
  calculated.}

  \item{\dots}{Optional arguments passed to methods like
  'print'}
}
\value{
  The summary function returns a list with the same
  components as 'summary.mkinfit', and the additional
  components from the results of the optimization routine
  used. %%\item{diffs}{} %%\item{}{} %%\item{}{}
  %%\item{}{}
}
\description{
  S3 method for class 'kingui0'
}
\note{
  This is deprecated!
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
summary(fit1)
fit2 <- IRLSkinfit.gui(SFO_SFO_gui)
summary(fit2)
}
}
\author{
  Zhenglei Gao
}
\seealso{
  \code{\link{summary.mkinfit}},
  \code{\link{summary.modFit}}
}
\keyword{internal}

