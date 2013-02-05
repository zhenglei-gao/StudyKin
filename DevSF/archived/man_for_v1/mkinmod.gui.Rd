\name{mkinmod.gui}
\alias{mkinmod.gui}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to set up a kinetic model with one or more compartments.
}
\description{
GUI version of \code{\link{mkinmod}}. The function takes a
specification, consisting of a list of the compartments in the
data. Each compartment is again represented by a list, specifying the
kinetic model type, reaction or transfer to other observed compartments,
the initial parameter values, lower and upper bounds, fixed or not, and
observed data.
}
\usage{
mkinmod.gui(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    Each list cell represents a comparment which contains a list of
    comonents including 'type'(kinetic reaction type, single first order
    kinetics "SFO" are implemented for all compartments, while "FOMC",
    "DFOP" and "HS" can additionally be chosen for the first variable
    which is assumed to be the source compartment),  each parameter
    name(a list of 'ini','fixed','lower','upper'),'residue'(measured
    concentrations),'time'(sampling time),'weight'(weights to be used,
    defaul 1),'sink'( Default TRUE, tranformation to unspecified
    compartments.),'to'(a vector of compartment names that the source
    compartment will be transferred to).
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
A list of class 'mkinmod.gui' for use with
\code{\link{mkinfit.gui}},\code{\link{IRLSkinfit.gui}} and
\code{\link{mcmckinfit.gui}} containing:

\item{diffs}{ A vector of string representations of differential equations,
          one for each modelling compartment.}
\item{parms}{ A vector of parameter names occurring in the differential
          equations.}
\item{map}{A list containing named character vectors for each
          compartments in the model.}
\item{parms.ini}{Initial values for all kinetic parameters in the
          model.}
\item{state.ini}{Initial state values for all compartments in the
          model.}
\item{lower}{Lower bounds for the parameters(including state variables)
          to be optimized.}
\item{upper}{upper bounds for the parameters(including state variables)
          to be optimized.}
\item{fixed_parms}{The names of the kinetic parameters that are fixed
          during optimization.}
\item{fixed_initials}{ The names of the initial states that are fixed
          during optimization.}
\item{residue}{The observed data matrix with a time column.}
\item{weightmat}{The weights matrix.}
\item{ff}{A vector of string representations of the transformation
between the formation fractions in the model and the transfomed
formation fractions in the optimization process.}
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
  \code{\link{mkinmod}}
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
12.7, 11.4, 11.6, 10.9,  9.5,  7.6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model }
% \keyword{  }% __ONLY ONE__ keyword per line
