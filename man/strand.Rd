\name{strand}
\docType{methods}
\alias{strand-methods}

\alias{strand}
\alias{strand,missing-method}
\alias{strand,character-method}

\title{Accessing strand information}

\description{The \code{strand} generic is meant as an accessor for
  strand information. Two methods are defined by the \code{BSgenome}
  package, described below.
}

\usage{strand(x, ...)}

\arguments{
  \item{x}{The object from which to obtain a strand factor, can be
    missing.}
  \item{...}{Additional arguments to pass to methods}
}

\details{
  If \code{x} is missing, returns an empty factor with the
  standard levels that any strand factor should have: \code{+},
  \code{-}, and \code{*} (for either).

  If \code{x} is a \code{character} vector, \code{x} is coerced to a
  factor with the levels listed above.
}

\author{ Michael Lawrence }

\examples{
strand()
strand(c("+", "-", NA, "*"))
}

\keyword{methods}