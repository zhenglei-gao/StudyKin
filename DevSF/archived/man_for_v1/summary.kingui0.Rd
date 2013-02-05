\name{summary.kingui}
\alias{summary.kingui}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 method for class 'kingui'
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{
summary.kingui(object, data = TRUE, distimes = TRUE, ff = TRUE, cov =
FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    %%     ~~Describe \code{object} here~~
    An object of class 'kingui' from the result of NLS or IRLS fit.
}
  \item{data}{
    %%     ~~Describe \code{data} here~~
    If TRUE, include in the returned values a data frame containing the
    observed and predicted values with residuals and estimated standard
    deviations or weights.
}
  \item{distimes}{
    %%     ~~Describe \code{distimes} here~~
    If TRUE, DT50 and DT90 values should be
          included.
}
  \item{ff}{
    %%     ~~Describe \code{ff} here~~
    If TRUE, the formation fraction should be calculated from the
    estimated transformed parameters.
}
  \item{cov}{
    %%     ~~Describe \code{cov} here~~
    If TRUE, parameter covariances should be
          calculated.
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    Optional arguments passed to methods like 'print'
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The summary function returns a list with the same components as
'summary.mkinfit', and the additional components from the results of the
optimization routine used.
%%\item{diffs}{}
%%\item{}{}
%%\item{}{}
%%\item{}{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{summary.mkinfit}}, \code{\link{summary.modFit}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Summary Statistics and Plots }
