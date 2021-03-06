\name{kinplot}
\alias{kinplot}
\title{Plot the fit results}
\usage{
  kinplot(fit, name = NULL, xlab = "Time",
    ylab = "Observed",
    xlim = c(1, 1.1) * range(fit$data$time),
    ylim = c(1, 1.1) * range(fit$data$observed, na.rm = TRUE),
    legend = TRUE, pdf = FALSE, ...)
}
\arguments{
  \item{fit}{An object of class 'kingui'.}

  \item{name}{Names of the compartmens that should be
  plotted.}

  \item{xlab}{xlab name.}

  \item{ylab}{ylab name.}

  \item{xlim}{The plot range for x axis.}

  \item{ylim}{The plot range for y axis.}

  \item{legend}{whether to include legends.}

  \item{pdf}{whether to put the plots into one pdf file.}

  \item{...}{optional arguments. Not used for now.}
}
\value{
  \code{NULL}
}
\description{
  Function to plot the fit results for the whole model or
  selected compartments.
}
\note{
  This function makes the plots.
}
\author{
  Zhenglei Gao
}
\seealso{
  \code{\link{kingraph}}
}

