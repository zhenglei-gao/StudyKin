\name{kingui.control}
\alias{kingui.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Choose the algorithm to use and related control parameters in kinetic
  evaluations using NLS and IRLS.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Allow the user to set some characteristics of the optimization procedure.
}
\usage{
kingui.control(method = "solnp", maxIter = 100, tolerance = 1e-08,
odesolver = "lsoda", atol = 1e-09, rtol = 1e-10, rhobeg = 0.05, iprint =
1, trace = 0, goMarq = 0, delta = 1e-06, rho = 1, submethod = "Marq",
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    %%     ~~Describe \code{method} here~~
    The method to be used, one of "Rsolnp", "minqa","spg", "Marq", "Port", "Newton",
          "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Pseudo" -
          see details.
}
  \item{maxIter}{
    %%     ~~Describe \code{maxIter} here~~
    Numerber of maximum iterations
}
  \item{tolerance}{
    %%     ~~Describe \code{tolerance} here~~
    A positive numeric value specifying the tolerance level for
          the relative offset convergence criterion.
}
  \item{odesolver}{
    %%     ~~Describe \code{odesolver} here~~
    the integration routines used. see \code{\link{ode}}
}
  \item{atol}{absolute error tolerance, either a scalar or an array as long
          as 'y'. See \code{\link{ode}}
%%     ~~Describe \code{atol} here~~
}
  \item{rtol}{
    %%     ~~Describe \code{rtol} here~~
    relative error tolerance, either a scalar or an array as long
          as 'y'.See \code{\link{ode}}.
}
  \item{rhobeg}{
    %%     ~~Describe \code{rhobeg} here~~
    Initial trust region radius for method 'bobyqa'. See details of see \code{\link{bobyqa}}.
}
  \item{iprint}{
    %%     ~~Describe \code{iprint} here~~
    control the amount of printing by setting iprint to 0, 1 2, or 3.  See details of see \code{\link{bobyqa}}.
}
  \item{trace}{
    %%     ~~Describe \code{trace} here~~
    A logical variable (TRUE/FALSE). If 'TRUE', information on
          the progress of optimization is printed. See details of see \code{\link{spg}}.
}
  \item{goMarq}{
    %%     ~~Describe \code{goMarq} here~~
    If TRUE, using "Marq" for the iterations after the first step in IRLS.
}
  \item{delta}{
    %%     ~~Describe \code{delta} here~~
    Control parameters for 'solnp'. See details of see \code{\link{solnp}}.
}
  \item{rho}{
    %%     ~~Describe \code{rho} here~~
    Control parameters for 'solnp'. See details of see \code{\link{solnp}}.
}
  \item{submethod}{
    %%     ~~Describe \code{submethod} here~~
 If the method chosen failed to produce results, run the optimization
  using a substitute method.
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    Other characteristics of different optimizer for the users to play
  with.
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
A list of control parameters for the ODE solver and the optimization
routines.
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
kingui.control()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Kinetic Evaluations}
\keyword{Optimization}% __ONLY ONE__ keyword per line
