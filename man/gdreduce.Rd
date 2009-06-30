\name{gdreduce}
\Rdversion{1.1}
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
gdreduce(f, ..., init, right = FALSE, accumulate = FALSE, gdArgs = list()) 
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
	dicteated by methods defined on \code{gdreduce}
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
\value{

  An object of class \code{GenomeData}, containing elements
  corresponding to the intersection of all named elements of
  \code{\dots}.

}
\seealso{
  \code{Reduce}
}
\author{
  Martin Morgan
}

\examples{
showMethods(gdreduce, where=getNamespace("BSgenome"))
}
\keyword{manip}

