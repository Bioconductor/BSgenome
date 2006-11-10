\name{getSeq}

\alias{getSeq}

\title{getSeq}

\description{
   Extract the sequence delimitted by start and end
   in a chromosome of a given genome.
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

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link{BSgenome-class},
  \code{\link{seqnames}},
  \code{\link[Biostrings]{BStringViews}}
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
}

\keyword{manip}

