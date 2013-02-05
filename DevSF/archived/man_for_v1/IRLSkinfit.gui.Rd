\name{IRLSkinfit.gui}
\alias{IRLSkinfit.gui}
\alias{summary.kingui}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a kinetic model using the IRLS algorithm.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
   This function does kinetic evaluations using the IRLS algorithm.
}
\usage{
IRLSkinfit.gui(mkinmodini, eigen = FALSE, plot = FALSE, plottitle = "",
quiet = FALSE, err = NULL, weight = "none", scaleVar = FALSE,
 ctr = kingui.control(), irls.control = list(), update = NULL, ...)

\method{summary}{kingui}(object, data = TRUE, distimes = TRUE, ff =
TRUE, cov = FALSE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mkinmodini}{
    %%     ~~Describe \code{mkinmodini} here~~
     A list of class \code{\link{mkinmod.gui}}, containing the kinetic
     model to be fitted to the data, and the initial parameter values,
     the observed data.
}
  \item{eigen}{
    %%     ~~Describe \code{eigen} here~~
    If TRUE,  the solution of the system of differential equations
    should be based on the spectral decomposition of the coefficient
    matrix in cases that this is possible.
}
  \item{plot}{
    %%     ~~Describe \code{plot} here~~
    If TRUE,the observed values and the numerical solutions should be plotted at each stage
    of the optimisation.
}
  \item{plottitle}{
    %%     ~~Describe \code{plottitle} here~~
    The title of the plot for visualizing the optimization process.
}
  \item{quiet}{
    %%     ~~Describe \code{quiet} here~~
     If TRUE, suppress printing out the current model cost after each(>1) improvement.
}
  \item{err}{
    %%     ~~Describe \code{err} here~~
    See argumetns of \code{\link{mkinfit.gui}}
}
  \item{weight}{
    %%     ~~Describe \code{weight} here~~
    See argumetns of \code{\link{mkinfit.gui}}
}
  \item{scaleVar}{
    %%     ~~Describe \code{scaleVar} here~~
    See argumetns of \code{\link{mkinfit.gui}}
}
  \item{ctr}{
    %%     ~~Describe \code{ctr} here~~
      a list of control values for the estimation algorithm to replace
     the default values including maximum iterations and absolute error
     tolerance.  Defaults to the output of \code{\link{kingui.control}}.
}
  \item{irls.control}{
    %%     ~~Describe \code{irls.control} here~~
    A list of control values for the estimation algorithm to replace
     the default values including the maximum number of iterations for
     the outer iteration and the error tolerance level for the error
     variance estimation updating.
}
  \item{update}{
    %%     ~~Describe \code{update} here~~
    If not NULL, should be a list of starting values obtained from other
    optimization methods.
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    Further arguments that will be passed to \code{\link{modFit}}.
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
A list with  "kingui", "mkinfit" and "modFit" in the class
  attribute. A summary can be obtained by \code{\link{summary.kingui}}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  Zhenglei Gao
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
   \code{\link{IRLSkinfit}},  \code{\link{mkinfit.gui}}
}
\examples{
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
fit <- IRLSkinfit.gui(SFO_SFO_gui,plot=T,ctr=kingui.control(maxIter=100,
tolerance=1e-6,odesolver='lsoda'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IRLS }
\keyword{ Kinetic Evaluations }% __ONLY ONE__ keyword per line
