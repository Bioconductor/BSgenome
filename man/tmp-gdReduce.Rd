\name{gdReduce}
\Rdversion{1.1}
\alias{gdReduce}
\alias{gdReduce,GenomeData-method}
\alias{gdReduce,GenomeDataList-method}

% Deprecated
\alias{gdreduce}

\title{
  Reduces arguments to a single GenomeData instance
}
\description{
  This function accepts one or more objects that are reduced, with a
  user-specified function, to a single \code{\linkS4class{GenomeData}}
  instance.

}
\usage{
gdReduce(f, ..., init, right = FALSE, accumulate = FALSE, gdArgs = list()) 
}
\arguments{
  \item{f}{
    An object of class \code{"function"}, accepting two instances of
	classes appropriate for the \code{...} arguments, and returning an
	object suitable for subsequent use in \code{f} and incorporation
	into \code{GenomeData}.
  }
  \item{\dots}{
	Objects to be reduced. All objects should be of the same class, as
	dictated by methods defined on \code{gdReduce}
    A function to be applied to each chromosome-level sub-element of
    \code{X}.
  }
  \item{init}{
	An R object of the same kind as the elements of \dots.
  }
  \item{right}{
	A logical indicating whether to proceed from left to right (default)
	or right to left.
  }
  \item{accumulate}{
	A logical indicating whether the successive reduce combinations
	should be accumulated. By default, only the final combination is
	used.}
  \item{gdArgs}{
	Additional arguments passed to the \code{\linkS4class{GenomeData}}
	constructor used to assemble the final object.
  }
}
\details{
  The \code{gdReduce} method for \link{GenomeData} objects successively
  combines \link{GenomeData} elements of \code{...} using \code{f};
  all arguments assigned to \code{...} must be of class \link{GenomeData}.
  \code{f} is a function accepting two objects returned by \code{"[["}
  applied to the successive elements of \code{...}, returning a single
  \link{GenomeData} object to be used in subsequent calls to \code{f}.
  \code{init}, \code{right}, and \code{accumulate} are as described for
  \code{\link{Reduce}}. \code{gdArgs} can be used to provide metadata
  information to the constructor used to create the final \link{GenomeData}
  object.

  Currently the \code{gdReduce} method for \link{GenomeDataList} objects
  works when a single \link{GenomeDataList} object \code{x} is provided
  as \code{...} and it does \code{gdReduce(f, x[[1]], x[[2]] ... x[[N]],
  init, right, accumulate, gdArgs)} where \code{N} is the length of \code{x}
  i.e. the number of \link{GenomeData} objects in it.
}
\value{
  An object of class \code{GenomeData}, containing elements
  corresponding to the intersection of all named elements of
  \code{\dots} (in the case of the method for \link{GenomeData} objects)
  or all elements in the single \link{GenomeDataList} object passed to
  it (in the case of the method for \link{GenomeDataList} objects).
}
\author{
  Martin Morgan
}
\seealso{
  \code{\link{Reduce}},
  \link{GenomeData-class},
  \link{GenomeDataList-class}
}
\examples{
gdReduce
showMethods("gdReduce")

gd <- GenomeData(list(chr1 = IRanges(1, 10), chrX = IRanges(2, 5)),
                 organism = "Mmusculus", provider = "UCSC",
                 providerVersion = "mm9")

gdr <- gdReduce(function(x, y) {
    ## "[[" returns IRanges instances, construct a synthetic version
    IRanges(c(start(x), start(y)), c(end(x), end(y)))
}, GenomeDataList(list(gd, gd[2])))
gdr[["chr1"]]
gdr[["chrX"]]
}
\keyword{manip}
