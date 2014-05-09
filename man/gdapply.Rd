\name{gdapply}
\Rdversion{1.1}
\alias{gdapply}
\alias{gdapply,GenomeData,function-method}
\alias{gdapply,GenomeDataList,function-method}

\title{
  Applies a function to elements of a GenomeData
}
\description{
  WARNING: Starting with BioC 3.0, GenomeData and GenomeDataList objects
  are deprecated. Note that the GenomeData/GenomeDataList containers
  predate the \link[GenomicRanges]{GRanges}/\link[GenomicRanges]{GRangesList}
  containers and, most of the times, the latters can be used instead of
  the formers.
  Please let us know on the bioc-devel mailing list
  (\url{http://bioconductor.org/help/mailing-list/}) if you have
  a use case where you think there are significant benefits in
  using GenomeData/GenomeDataList over
  \link[GenomicRanges]{GRanges}/\link[GenomicRanges]{GRangesList},
  or if you have questions or concerns about this.

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
