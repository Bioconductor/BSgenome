\name{getSeq}

\alias{getSeq}

\title{getSeq}

\description{
   Extract the sequence delimitted by start and end
   in a chromosome of a given genome.
}

\usage{
  getSeq(bsgenome, seqname, start=NA, end=NA, as.character=TRUE)
}

\arguments{
  \item{bsgenome}{
  }
  \item{seqname}{
  }
  \item{start}{
  }
  \item{end}{
  }
  \item{as.character}{
  }
}

\details{
}

\value{
  By default, a character vector.
  If \code{as.character=FALSE}, a \code{\link[Biostrings]{BStringViews}} object.
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class},
  \code{\link{seqnames}},
  \code{\link[Biostrings]{BStringViews}}
}

\examples{
  # Load the Caenorhabditis elegans genome (UCSC Release ce2):
  library(BSgenome.Celegans.UCSC.ce2)

  # Look at the index of sequences:
  Celegans

  # Get 10-base sequences starting at position 150, 250, ..., 2050
  # in chromosome V:
  starts <- seq(150, by=100, length.out=20)
  getSeq(Celegans, "chrV", starts, starts + 9)

  # The same returned as a BStringViews object (_much_ faster, no data is copied):
  getSeq(Celegans, "chrV", starts, starts + 9, as.character=FALSE)

  # The display of a BStringViews object stays compact:
  starts <- seq(-10, by=1000, length.out=50)
  getSeq(Celegans, "chrV", starts, starts + 1002, as.character=FALSE)
  # Note that the views overlap.
}

\keyword{manip}

