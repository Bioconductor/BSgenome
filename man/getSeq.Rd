\name{getSeq}

\alias{getSeq}
\alias{getSeq,BSgenome-method}

\title{getSeq}

\description{
   A convenience function for extracting a set of sequences (or
   subsequences) from a \link{BSgenome} or other object. This man page
   specifically documents the \link{BSgenome} method.
}

\usage{
  \S4method{getSeq}{BSgenome}(bsgenome, names, start=NA, end=NA, width=NA, strand="+", as.character=TRUE)
}

\arguments{
  \item{bsgenome}{
    A \link{BSgenome} object.
    See the \code{\link{available.genomes}} function for how to install
    a genome.
  }
  \item{names}{
    The names of the sequences to extract from \code{bsgenome}.
    If missing, then \code{seqnames(bsgenome)} is used.

    See \code{?\link{seqnames}} and \code{?\link{mseqnames}} to get
    the list of single sequences and multiple sequences (respectively)
    contained in \code{bsgenome}.

    Here is how the lookup between the names passed to the \code{names}
    argument and the sequences in \code{bsgenome} is performed.
    For each \code{name} in \code{names}:
    (1) if \code{bsgenome} contains a single sequence with that name
    then this sequence is returned;
    (2) otherwise the names of all the elements in all the multiple
    sequences are searched: \code{name} is treated as a regular
    expression and \code{\link[base]{grep}} is used for this search.
    If exactly one sequence is found, then it's returned, otherwise an
    error is raised.
  }
  \item{start, end, width}{
    Vector of integers (eventually with NAs).
  }
  \item{strand}{
    A vector containing \code{+}s or/and \code{-}s.
  }
  \item{as.character}{
    \code{TRUE} or \code{FALSE}. Should the extracted sequences
    be returned in a standard character vector?
  }
}

\details{
  The \code{names}, \code{start}, \code{end}, \code{width} and \code{strand}
  arguments are expanded cyclically to the length of the longest provided
  none are of zero length.
}

\value{
  A standard character vector when \code{as.character=TRUE}.
  Note that when \code{as.character=TRUE}, then the masks that
  are defined on top of the sequences to extract are ignored
  (i.e. dropped) if any (see \code{?`\link[Biostrings]{MaskedXString-class}`}
  for more information about masked sequences).

  A \link[Biostrings]{DNAString} or \link[Biostrings]{MaskedDNAString}
  object when \code{as.character=FALSE}.
  Note that \code{as.character=FALSE} is not supported yet when extracting
  more than one sequence.
}

\note{
  Be aware that using \code{as.character=TRUE} can be very inefficient
  when the returned character vector contains very long strings
  (> 1 million letters) or is itself a long vector (> 10000 strings).

  \code{getSeq} is much more efficient when used with
  \code{as.character=FALSE} but this works only for extracting
  one sequence at a time for now.
}

\author{H. Pages; improvements suggested by Matt Settles}

\seealso{
  \code{\link{available.genomes}},
  \link{BSgenome-class},
  \code{\link{seqnames}},
  \code{\link{mseqnames}},
  \code{\link[base]{grep}},
  \code{\link[IRanges]{subseq}},
  \code{\link[Biostrings]{DNAString}},
  \code{\link[Biostrings]{MaskedDNAString}},
  \code{\link{[[,BSgenome-method}}
}

\examples{
  # Load the Caenorhabditis elegans genome (UCSC Release ce2):
  library(BSgenome.Celegans.UCSC.ce2)

  # Look at the index of sequences:
  Celegans

  # Get chromosome V as a DNAString object:
  getSeq(Celegans, "chrV", as.character=FALSE)
  # which is in fact the same as doing:
  Celegans$chrV

  # Never try this:
  #getSeq(Celegans, "chrV")
  # or this (even worse):
  #getSeq(Celegans)

  # Get the first 20 bases of each chromosome:
  getSeq(Celegans, end=20)

  # Get the last 20 bases of each chromosome:
  getSeq(Celegans, start=-20)

  # Extracting small sequences from different chromosomes:
  myseqs <- data.frame(
    chr=c("chrI", "chrX", "chrM", "chrM", "chrX", "chrI", "chrM", "chrI"),
    start=c(NA, -40, 8510, 301, 30001, 9220500, -2804, -30),
    end=c(50, NA, 8522, 324, 30011, 9220555, -2801, -11),
    strand=c("+", "-", "+", "+", "-", "-", "+", "-")
  )
  getSeq(Celegans, myseqs$chr,
         start=myseqs$start, end=myseqs$end)
  getSeq(Celegans, myseqs$chr,
         start=myseqs$start, end=myseqs$end, strand=myseqs$strand)

  # Get the "NM_058280_up_1000" sequence (belongs to the upstream1000
  # multiple sequence) as a character string:
  s1 <- getSeq(Celegans, "NM_058280_up_1000")
  # or a DNAString object (more efficient):
  s2 <- getSeq(Celegans, "NM_058280_up_1000", as.character=FALSE)

  getSeq(Celegans, "NM_058280_up_5000", start=-1000) == s1  # TRUE

  getSeq(Celegans, "NM_058280_up_5000",
         start=-1000, as.character=FALSE) == s2  # TRUE
}

\keyword{manip}
