\name{completeCompound}
\alias{completeCompound}
\title{Auxiliary function including error checks for \code{mkinmod.full}}
\usage{
  completeCompound(compound = list(type = "SFO", to = "M1"),
    varname = NULL, first = FALSE,
    inpartri = c("default", "water-sediment", "advanced"),
    outpartri = c("default", "water-sediment", "advanced"),
    data = NULL, weight = NULL, update = NULL, ...)
}
\arguments{
  \item{compound}{A list of properties for a single
  compound.}

  \item{varname}{The compound name.}

  \item{first}{hether this compound is the parent compound}

  \item{inpartri}{Input parameterization.}

  \item{outpartri}{Output parameterization.}

  \item{data}{If not NULL, The residue data frame.}

  \item{weight}{If weight is NULL, check weight component.}

  \item{update}{If not NULL, replace the components in
  compound with the ones in the update.}

  \item{...}{Other optional arguments. Not used.}
}
\value{
  A list of components to be used in
  \code{\link{mkinmod.full}}. For the differential
  functions related to other compounds, they are derived in
  the \code{\link{mkinmod.full}} function.
}
\description{
  Auxiliary function including error checks for
  \code{mkinmod.full}
}
\author{
  Zhenglei Gao
}

