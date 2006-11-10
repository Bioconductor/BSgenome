\name{getSeq}

\alias{getSeq}

\title{getSeq}

\description{
   Extract the sequence delimitted by a start and an end positions
   relative to the full sequence of a given genome.
}

\usage{
  getSeq(bsgenome, seqname, start=NA, end=NA, as.BStringViews=FALSE)
}

\arguments{
  \item{bsgenome}{
    A \link{BSgenome} object as one of those found in the BSgenome
    data packages. See the \code{\link{available.genomes}} function
    for more details.
  }
  \item{seqname}{
    [TODO: Document me]
  }
  \item{start}{
    [TODO: Document me]
  }
  \item{end}{
    [TODO: Document me]
  }
  \item{as.BStringViews}{
    [TODO: Document me]
  }
}

\details{
}

\value{
  If \code{length(seqname) == 1}, a character vector (\code{as.BStringViews=FALSE})
  or a \code{\link[Biostrings]{BStringViews}} object (\code{as.BStringViews=TRUE}).
  If \code{length(seqname) != 1}, a character vector (\code{as.BStringViews}
  is ignored).
}

\note{
  \code{getSeq} is very efficient when used with \code{as.BStringViews=TRUE}
  (this only works when \code{seqname} is of length 1), because in this case
  the sequence data are not copied.
  Otherwise, the data are copied. Be aware that this can be very inefficient
  if the returned vector contains very long strings (> 1 million letters)
  or is itself very long (> 10000 strings).
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link{BSgenome-class},
  \code{\link{seqnames}},
  \code{\link[base]{substr}},
  \code{\link[Biostrings]{subBString}},
  \code{\link[Biostrings]{BStringViews}},
  \code{\link[Biostrings]{DNAString}}
}

\examples{
  # Load the Caenorhabditis elegans genome (UCSC Release ce2):
  library(BSgenome.Celegans.UCSC.ce2)

  # Look at the index of sequences:
  Celegans

  # Get the first 20 bases of each chromosome:
  getSeq(Celegans, seqnames(Celegans), 1, 20)

  # Some sequences starting at pos 1 in chromosome V:
  getSeq(Celegans, "chrV", 1, 1:4)

  # Omitting the 'start' (or the 'end') argument is equivalent
  # to starting at the first nucleotide (or ending at the last
  # nucleotide):
  getSeq(Celegans, "chrV", , 1:4)
  getSeq(Celegans, "chrV", 20922200)

  # Never try this:
  #getSeq(Celegans, "chrV")
  # or this (even worse):
  #getSeq(Celegans, seqnames(Celegans))
  # unless you want to see millions of screens filled with A, C, G and T
  # and kill your system.

  # Get the 10-base sequences starting at positions 150, 250, ..., 1250
  # in chromosome V:
  starts <- seq(150, by=100, length.out=12)
  getSeq(Celegans, "chrV", starts, starts + 9)

  # The same returned as a BStringViews object (_much_ faster, no data is copied):
  getSeq(Celegans, "chrV", starts, starts + 9, as.BStringViews=TRUE)

  # The display of a BStringViews object stays compact:
  starts <- seq(-10, by=1000, length.out=50)
  getSeq(Celegans, "chrV", starts, starts + 1002, as.BStringViews=TRUE)
  # Note that those views have overlaps of 3 nucleotides.

  # Get the whole chromosome V sequence as a DNAString object:
  getSeq(Celegans, "chrV", as.BStringViews=TRUE)[[1]]
}

\keyword{manip}

