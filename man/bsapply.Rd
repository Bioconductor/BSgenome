\name{bsapply}

\alias{bsapply}

\title{bsapply}

\description{
  Apply a function to each chromosome in a genome.
}

\usage{
  bsapply(X, FUN, exclude, ...)
}

\arguments{
  \item{X}{
    a BSgenome object that contains chromosomes that you wish to apply
  FUN on
  }
  \item{FUN}{
    the function to apply to each chromosome in the BSgenome object
  }
  \item{exclude}{
    this is a character vector with strings that will be used to filter
  out chromosomes whose names match these strings.
  }
  \item{...}{
    optional arguments to 'FUN'.
  }
}

\details{
  By default the exclude parameter is set to not exclude anything.  A
  popular option will probably be to set this to "rand" so that random
  bits of unassigned contigs are filtered out.
}

\value{
  A named list is returned where each element is whatever was returned by the
  function FUN.  The names of the list correspond to the names of the
  chromosomes as they were labeled by the BSgenome obect in the slot
  accessed by seqnames().
}

\author{Marc Carlson}

\seealso{
  \link{BSgenome-class}
}

\examples{
  ## Load the Worm genome:
  library("BSgenome.Celegans.UCSC.ce2")

  ## Count the alphabet frequencies for every chromosome but exclude
  ## mitochrondrial ones:
  bsapply(Celegans, alphabetFrequency, exclude = "M")

}

\keyword{manip}
