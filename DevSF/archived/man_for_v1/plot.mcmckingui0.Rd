\name{plot.mcmckingui}
\alias{plot.mcmckingui}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  S3 method to plot for calss 'mcmckingui'
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.mcmckingui(object, fname1, fname2, pch = 1, device = "wmf", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    %%     ~~Describe \code{object} here~~
    An object of class 'mcmckingui'
}
  \item{fname1}{
    %%     ~~Describe \code{fname1} here~~
    The file name of the density plot.
}
  \item{fname2}{
    %%     ~~Describe \code{fname2} here~~
     The file name of the correlation plot.
}
  \item{pch}{
    %%     ~~Describe \code{pch} here~~
    What kind of points to use in the plots.
}
  \item{device}{
    %%     ~~Describe \code{device} here~~
    The plot device to be used.
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    Other arguments to be passed to 'plot'.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Density and Correlation plots of the sampled parameters in 'wmf' or other format.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
}
\examples{
\dontrun{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Summary Statistics and Plots }

