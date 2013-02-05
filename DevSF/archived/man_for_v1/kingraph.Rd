\name{kingraph}
\alias{kingraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate graph data for the GUI to make fit plot.
}
\description{
A data table is generated and stored in file=filename for later usage.
}
\usage{
kingraph(fit, filename = "graphdata.txt",xlim = c(1, 1.05) *
range(fit$data$time), ylim = c(1, 1.05)* range(fit$data$observed, na.rm
= TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    %%     ~~Describe \code{fit} here~~
    An object of class 'kingui'.
}
  \item{filename}{
    %%     ~~Describe \code{filename} here~~
    The file in which the graph data will be stored.
}

  \item{xlim}{
    %%     ~~Describe \code{xlim} here~~
    The plot range for x axis.
}
  \item{ylim}{
    %%     ~~Describe \code{ylim} here~~
    The plot range for y axis.
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
}
\examples{
\dontrun{
SFO_SFO_gui <- mkinmod.gui(Parent = list(type = "SFO", to = "Metab",
              sink = TRUE,
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
fit <- mkinfit.gui(SFO_SFO_gui,plot=T)
kingraph(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Summary Statistics and Plots }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
