\name{gdapply}
\Rdversion{1.1}
\alias{gdapply}
\alias{gdapply,GenomeData,function-method}
\alias{gdapply,GenomeDataList,function-method}

% Defunct:
\alias{gdApply}

\title{
  Applies a function to elements of a GenomeData
}
\description{
  Returns a list of values obtained by applying a function to elements
  of a \link{GenomeData} or \link{GenomeDataList} object.
}
\usage{
gdapply(X, FUN, ...)
}
\arguments{
  \item{X}{
    An object of class \link{GenomeData} or \link{GenomeDataList}.
  }
  \item{FUN}{
    A function to be applied to each chromosome-level sub-element of
    \code{X}.
  }
  \item{\dots}{
    Further arguments; passed to \code{FUN}
  }
}
\value{
  Typically an object of the same class as \code{X}.
}
\author{
  Deepayan Sarkar
}
\seealso{
  \link{GenomeData-class},
  \link{GenomeDataList-class}
}
\keyword{manip}
